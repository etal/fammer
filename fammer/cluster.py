"""Extract phylogenetically-based clusters from a sequence set.

Uses FastTree to quickly build a tree with branch support values, then extracts
well-supported clades from the tree and writes the corresponding sequence sets
to FASTA files. Unclustered sequences are written to another "Unique" file.
Also writes a phyloXML tree file (.xml) showing clusters as colorized clades.
The output files are written to the current directory and named after the
base name of the input file.

Depends:
    biopython (pip install biopython)
    FastTree, as fasttree or FastTree (apt-get install fasttree)
    mafft (apt-get install mafft)

"""

import colorsys
import logging
import math
import os.path
import subprocess
import tempfile
from collections import Counter
from io import StringIO

from Bio import AlignIO, Phylo, SeqIO
from Bio.Phylo.BaseTree import BranchColor
from Bio.Seq import Seq
from biofrills import alnutils

from ._share import UsageError, ungapped, write_fasta
from .tasks import which

logger = logging.getLogger(__name__)

SUFFIX_ERROR = (
    "Input sequences must be a Clustal alignment (.aln) or unaligned FASTA (.fasta)"
)

# Debian installs FastTree as fasttree, and bioconda as FastTree.
FASTTREE_NAMES = ("fasttree", "FastTree")


def cmd_cluster(args):
    """Cut on longest branches, then group by remaining supported clusters."""
    # Args:
    #     sequences: Target sequence file (foo.aln, clustal format)
    # Opts (ENH):
    #     --size: MIN_CLUSTER_SIZE (e.g. 4)
    #     --confidence: MIN_CONFIDENCE (e.g. 0.7, 0.9) -- needed?
    #     --tree: precomputed tree, incl. supports -- good idea?
    #     --outgroup: needed?
    if not os.path.isfile(args.sequences):
        raise UsageError(f"No sequence file at {args.sequences}")
    if not args.sequences.endswith((".aln", ".fasta")):
        raise UsageError(f"{SUFFIX_ERROR}: {args.sequences}")
    seq_idx = read_sequences(args.sequences)
    tree = load_tree(args.sequences)

    # Extract clusters of acceptable size
    MIN_CLUSTER_SIZE = max(4, int(round(tree.count_terminals() ** 0.4)))
    clusters = {}
    unclustered = set()
    seen_subnodes = set([tree.root])

    # Cluster according to long internal branches first
    # Take clades below "cut" branches as clusters
    for clade in get_clades_to_cut(tree):
        subnodes = list(clade.find_clades())
        if not seen_subnodes.isdisjoint(subnodes):
            continue
        taxa = set([k.name for k in subnodes if k.is_terminal()])
        if len(taxa) >= MIN_CLUSTER_SIZE:
            clusters[clade] = taxa
            seen_subnodes.update(subnodes)

    # Find any remaining supported clades & assign those to clusters too
    # Any tips not assigned to a cluster -> Unique
    for clade in tree.find_clades(order="level"):
        if clade in seen_subnodes:
            continue
        elif clade.is_terminal():
            unclustered.add(clade.name)
        else:
            subnodes = list(clade.find_clades())
            # If any part of this clade was a "cut" cluster, skip down to the
            # next layer to find untouched sub-clades
            if not seen_subnodes.isdisjoint(subnodes):
                continue
            # Extract all terminals into a new cluster
            seen_subnodes.update(subnodes)
            taxa = set([k.name for k in subnodes if k.is_terminal()])
            if len(taxa) >= MIN_CLUSTER_SIZE:
                clusters[clade] = taxa
            else:
                unclustered.update(taxa)

    write_clusters(args.sequences, seq_idx, tree, clusters, unclustered)


def read_sequences(seqfname):
    """The input records by identifier, ungapped if read from an alignment.

    Raises UsageError if an identifier occurs more than once, before MAFFT or
    FastTree is run, since the clusters could not name such a sequence.
    """
    is_aln = seqfname.endswith(".aln")
    with SeqIO.parse(seqfname, "clustal" if is_aln else "fasta") as parsed:
        records = list(parsed)
    counts = Counter(rec.id for rec in records)
    duplicates = sorted(seqid for seqid, count in counts.items() if count > 1)
    if duplicates:
        raise UsageError(
            f"Duplicate sequence identifiers in {seqfname}: " + ", ".join(duplicates)
        )
    if is_aln:
        for rec in records:
            rec.seq = Seq(ungapped(rec.seq))
    return {rec.id: rec for rec in records}


def load_tree(seqfname):
    """Load an alignment, build & prep a tree, return the tree object."""
    if not seqfname.endswith((".aln", ".fasta")):
        # A built-in ValueError, not UsageError: tests/data/py27/library/
        # cluster.json records this type, and the reference comparison looks
        # exception types up among the built-in names.
        raise ValueError(SUFFIX_ERROR)
    fasttree = fasttree_program()
    if seqfname.endswith(".aln"):
        aln = AlignIO.read(seqfname, "clustal")
    else:
        # Run MAFFT quickly
        alndata = subprocess.check_output(
            ["mafft", "--quiet", "--auto", seqfname], text=True
        )
        aln = AlignIO.read(StringIO(alndata), "fasta")

    # Use conserved (less-gappy) blocks to build the tree
    aln = alnutils.blocks(aln, 0.4)
    # FastTree ends a name at ':', ',', a parenthesis or whitespace, so the
    # identifiers it writes back can be truncated and no longer unique.  It
    # is given placeholder names instead, and the tips of the tree it returns
    # are renamed to the original identifiers.
    seqids = {}
    for index, rec in enumerate(aln):
        name = f"seq{index}"
        seqids[name] = rec.id
        rec.id = name
        rec.description = ""
    with tempfile.NamedTemporaryFile(mode="w") as tmp:
        AlignIO.write(aln, tmp, "fasta")
        tmp.flush()
        treedata = subprocess.check_output(
            [fasttree, "-pseudo", "-gamma", "-wag", tmp.name], text=True
        )
    tree = Phylo.read(StringIO(treedata), "newick")
    for tip in tree.get_terminals():
        tip.name = seqids[tip.name]
    return prepare_tree(tree)


def fasttree_program():
    """The path of the FastTree executable, under either of its names."""
    path = which(FASTTREE_NAMES)
    if path is None:
        raise UsageError(
            "FastTree is not on PATH under either name: " + ", ".join(FASTTREE_NAMES)
        )
    return path


def prepare_tree(tree):
    """Collapse weakly supported splits, ladderize, and zero the root branch.

    A split is weakly supported when its support is below the mean of the
    tree's support values.  A tree without any support values gives no
    evidence that a split is weak, so nothing is collapsed.
    """
    confs = [c.confidence for c in tree.find_clades() if c.confidence is not None]
    if confs:
        # ENH: accept min_confidence as an option
        min_confidence = math.fsum(confs) / len(confs)
        # A clade without a support value counts as weakly supported.
        tree.collapse_all(
            lambda c: c.confidence is None or c.confidence < min_confidence
        )
    tree.ladderize(reverse=True)
    tree.root.branch_length = 0.0
    return tree


def get_clades_to_cut(tree):
    """Select long internal branches to cut the tree at.

    A branch is long when it is at least the geometric mean of the tip
    branch lengths.  A tip of length zero, which FastTree gives to a sequence
    identical to another, has no logarithm and is left out of the mean, and
    a branch of length zero is never long.  When every tip has length zero,
    every branch of positive length is long.
    """
    tip_lengths = [cl.branch_length for cl in tree.get_terminals() if cl.branch_length]
    tip_mean = geom_mean(tip_lengths) if tip_lengths else 0.0
    clades_to_cut = [
        cl
        for cl in tree.get_nonterminals()
        if cl.branch_length and cl.branch_length >= tip_mean
    ]
    clades_to_cut.sort(key=lambda cl: cl.branch_length, reverse=True)
    logger.info("# to cut (up to): %s", len(clades_to_cut))
    return clades_to_cut


def geom_mean(values):
    """Geometric mean of a non-empty sequence of positive numbers.

    It is computed in log space, because the product of a few hundred branch
    lengths underflows to zero.
    """
    return math.exp(math.fsum(map(math.log, values)) / len(values))


def write_clusters(seqfname, seq_idx, tree, clusters, unclustered):
    """Write output files: clusters & unique as FASTA, tree as phyloXML.

    seq_idx maps each identifier to its record, as read_sequences returns it.
    Clusters are numbered by decreasing size, and clusters of equal size by
    their sorted member identifiers.
    """

    def write_cluster(seqids, fname):
        """Write the sequences of the given identifiers to a FASTA file."""
        records = [seq_idx[seqid] for seqid in seqids]
        write_fasta(records, fname)
        logger.info("Wrote %s (%d sequences)", fname, len(records))

    ordered = sorted(
        ((sorted(taxa), clade) for clade, taxa in clusters.items()),
        key=lambda item: (-len(item[0]), item[0]),
    )
    for i, (seqids, clade) in enumerate(ordered):
        write_cluster(seqids, os.path.basename(seqfname) + "." + str(i))
        # Evenly spaced hues, one per cluster.
        rgb = colorsys.hsv_to_rgb(i / len(ordered), 0.7, 0.9)
        clade.color = BranchColor(*(int(x * 255) for x in rgb))
        clade.width = 2
    if unclustered:
        write_cluster(sorted(unclustered), os.path.basename(seqfname) + ".Unique")

    treefname = os.path.basename(seqfname) + ".xml"
    Phylo.write(tree, treefname, "phyloxml")
    logger.info("Wrote %s", treefname)
