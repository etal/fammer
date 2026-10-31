"""Extract phylogenetically-based clusters from a sequence set.

Uses FastTree to quickly build a tree with branch support values, or reads a
given Newick tree, then extracts well-supported clades from the tree and writes
the corresponding sequence sets to FASTA files. Unclustered sequences are
written to another "Unique" file.
Also writes a phyloXML tree file (.xml) showing clusters as colorized clades,
and logs the silhouette of the clusters, from the patristic distances of the
tree.
The output files are written to the current directory and named after the
base name of the input file.

Depends:
    biopython (pip install biopython)
    numpy (pip install numpy)
    FastTree, as fasttree or FastTree (apt-get install fasttree)
    mafft (apt-get install mafft)

"""

import colorsys
import logging
import math
import os.path
import subprocess
import tempfile
from io import StringIO

import numpy as np
from Bio import AlignIO, Phylo, SeqIO
from Bio.Phylo.BaseTree import BranchColor
from Bio.Phylo.NewickIO import NewickError
from Bio.Seq import Seq
from biofrills import alnutils

from ._share import UsageError, repeated_names, ungapped, write_fasta
from .tasks import which

logger = logging.getLogger(__name__)

SUFFIX_ERROR = (
    "Input sequences must be a Clustal alignment (.aln) or unaligned FASTA (.fasta)"
)

# Debian installs FastTree as fasttree, and bioconda as FastTree.
FASTTREE_NAMES = ("fasttree", "FastTree")


def cmd_cluster(args):
    """Cut on longest branches, then group by remaining supported clusters.

    args carries the sequence file name (sequences) and three options, each
    None for its default: the smallest cluster size (size), the support below
    which a split is collapsed (confidence), and a Newick tree to use instead
    of building one with MAFFT and FastTree (tree).
    """
    if not os.path.isfile(args.sequences):
        raise UsageError(f"No sequence file at {args.sequences}")
    if not args.sequences.endswith((".aln", ".fasta")):
        raise UsageError(f"{SUFFIX_ERROR}: {args.sequences}")
    seq_idx = read_sequences(args.sequences)
    if args.tree is None:
        tree = load_tree(args.sequences)
    else:
        tree = read_tree(args.tree, set(seq_idx))
    # The silhouette is measured on the tree as built.  Collapsing a split
    # adds its branch length to each child, which moves the tips below it
    # apart, so distances taken afterwards would depend on the confidence
    # threshold as well as the clusters do.
    names, distances = patristic_distances(tree)
    tree = prepare_tree(tree, args.confidence)
    clusters, unclustered = find_clusters(tree, args.size)

    ordered = ordered_clusters(clusters)
    basename = os.path.basename(args.sequences)
    write_clusters(basename, seq_idx, tree, ordered, unclustered)
    log_silhouettes(basename, ordered, unclustered, names, distances)


def read_sequences(seqfname):
    """The input records by identifier, ungapped if read from an alignment.

    Raises UsageError if an identifier occurs more than once, before MAFFT or
    FastTree is run, since the clusters could not name such a sequence.
    """
    is_aln = seqfname.endswith(".aln")
    with SeqIO.parse(seqfname, "clustal" if is_aln else "fasta") as parsed:
        records = list(parsed)
    duplicates = repeated_names(rec.id for rec in records)
    if duplicates:
        raise UsageError(
            f"Duplicate sequence identifiers in {seqfname}: " + ", ".join(duplicates)
        )
    if is_aln:
        for rec in records:
            rec.seq = Seq(ungapped(rec.seq))
    return {rec.id: rec for rec in records}


def load_tree(seqfname):
    """Align the sequences if need be, build a tree with FastTree, and return
    it with the tips named by the sequence identifiers."""
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
    return tree


def read_tree(treefname, seqids):
    """Read a Newick tree whose tips are named by exactly the given identifiers.

    A tip that names no sequence could not be written to a cluster, and a
    sequence that no tip names would be left out of every output file.  An
    identifier containing whitespace or any of the characters :,;()[]' must be
    quoted in the Newick file; unquoted, it is read as a shorter name.
    """
    if not os.path.isfile(treefname):
        raise UsageError(f"No tree file at {treefname}")
    try:
        tree = Phylo.read(treefname, "newick")
    except (ValueError, NewickError) as exc:
        raise UsageError(
            f"Could not read one Newick tree from {treefname}: {exc}"
        ) from exc
    names = [tip.name for tip in tree.get_terminals()]
    unnamed = sum(not name for name in names)
    if unnamed:
        raise UsageError(f"Tips of {treefname} without a name: {unnamed}")
    duplicates = repeated_names(names)
    if duplicates:
        raise UsageError(
            f"Duplicate tip names in {treefname}: " + ", ".join(duplicates)
        )
    present = set(names)
    unknown = sorted(present - seqids)
    if unknown:
        raise UsageError(
            f"Tips of {treefname} that name no input sequence: " + ", ".join(unknown)
        )
    missing = sorted(seqids - present)
    if missing:
        raise UsageError(
            f"Input sequences that no tip of {treefname} names: " + ", ".join(missing)
        )
    return tree


def fasttree_program():
    """The path of the FastTree executable, under either of its names."""
    path = which(FASTTREE_NAMES)
    if path is None:
        raise UsageError(
            "FastTree is not on PATH under either name: " + ", ".join(FASTTREE_NAMES)
        )
    return path


def prepare_tree(tree, min_confidence=None):
    """Collapse weakly supported splits, ladderize, and zero the root branch.

    The splits are the internal clades other than the root; collapse_all
    never collapses the root.  A split is weakly supported when its support
    is below min_confidence, or by default below the mean of the splits'
    support values; a split without a support value counts as weakly
    supported.  A tree whose splits have no support values gives no evidence
    that a split is weak, so by default nothing is collapsed, and an explicit
    min_confidence is an error.  So is one above every split's support, which
    would collapse the tree to a single node, and one above 0 and up to 1 on
    a tree whose supports run above 1, such as bootstrap percentages, which
    would collapse almost no split.  A threshold of 0, which collapses only
    the splits without a support value, means the same on either scale.
    """
    splits = [
        c.confidence
        for c in tree.get_nonterminals()
        if c is not tree.root and c.confidence is not None
    ]
    if min_confidence is None:
        if splits:
            min_confidence = math.fsum(splits) / len(splits)
    elif not splits:
        raise UsageError(
            "No split of the tree has a support value, so no confidence threshold "
            "applies"
        )
    elif min_confidence > max(splits):
        raise UsageError(
            f"No split of the tree has support of {min_confidence:g} or more; the "
            f"largest is {max(splits):g}, and FastTree writes support values "
            "between 0 and 1, while bootstrap percentages run from 0 to 100"
        )
    elif 0 < min_confidence <= 1 < max(splits):
        raise UsageError(
            f"A threshold of {min_confidence:g} would collapse almost no split of "
            f"this tree, whose support values run up to {max(splits):g}; give the "
            "threshold on the tree's own scale, such as 90 rather than 0.9 for "
            "bootstrap percentages"
        )
    if min_confidence is not None:
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


def find_clusters(tree, min_size=None):
    """The clusters of a prepared tree, and the tips left out of them.

    A cluster is a clade of at least min_size tips, by default the larger of
    4 and the number of tips to the power 0.4, rounded; the whole tree is
    never a cluster, and a tip clade is left out rather than taken as one.
    The clades below the long branches get_clades_to_cut selects are taken
    first, longest branch first; then each remaining clade that no cluster
    overlaps is taken in level order, as a cluster if it has enough tips, and
    otherwise its tips are left out.

    Returns the clusters, as a mapping of each clade to the names of its
    tips, and the set of names of the tips left out.
    """
    if min_size is None:
        min_size = max(4, round(tree.count_terminals() ** 0.4))
    clusters = {}
    unclustered = set()
    seen_subnodes = set([tree.root])
    for clade in get_clades_to_cut(tree):
        subnodes = list(clade.find_clades())
        if not seen_subnodes.isdisjoint(subnodes):
            continue
        taxa = set([k.name for k in subnodes if k.is_terminal()])
        if len(taxa) >= min_size:
            clusters[clade] = taxa
            seen_subnodes.update(subnodes)
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
            seen_subnodes.update(subnodes)
            taxa = set([k.name for k in subnodes if k.is_terminal()])
            if len(taxa) >= min_size:
                clusters[clade] = taxa
            else:
                unclustered.update(taxa)
    return clusters, unclustered


def ordered_clusters(clusters):
    """The clusters, a mapping of clades to member identifiers, as a list of
    (sorted member identifiers, clade) pairs in the order they are numbered:
    by decreasing size, and clusters of equal size by their members."""
    return sorted(
        ((sorted(taxa), clade) for clade, taxa in clusters.items()),
        key=lambda item: (-len(item[0]), item[0]),
    )


def patristic_distances(tree):
    """The tip names, and the matrix of patristic distances between the tips.

    The distance between two tips is the sum of their depths less twice the
    depth of their last common ancestor.  In post-order the tips below each
    clade are consecutive, so the pairs whose last common ancestor is a given
    clade are the blocks between the ranges of its children, and each pair is
    assigned once.
    """
    depths = tree.depths()
    count = tree.count_terminals()
    # The depth of the last common ancestor of each pair of tips.
    shared = np.zeros((count, count))
    names = []
    rooted = []
    spans = {}
    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            spans[clade] = (len(names), len(names) + 1)
            names.append(clade.name)
            rooted.append(depths[clade])
            continue
        children = [spans[child] for child in clade.clades]
        for index, (start, stop) in enumerate(children):
            for other_start, other_stop in children[index + 1 :]:
                shared[start:stop, other_start:other_stop] = depths[clade]
                shared[other_start:other_stop, start:stop] = depths[clade]
        spans[clade] = (children[0][0], children[-1][1])
    rooted = np.array(rooted)
    # The distances are written over the shared depths, to hold one matrix
    # fewer, as (depth_i + depth_j) - 2 * shared depth.
    shared *= 2.0
    distances = np.subtract(np.add.outer(rooted, rooted), shared, out=shared)
    np.fill_diagonal(distances, 0.0)
    return names, distances


def silhouette(distances, labels):
    """The silhouette of each point, from a square matrix of distances.

    labels gives the group of each point, and there must be at least two
    groups.  As in sklearn.metrics.silhouette_samples, a point alone in its
    group scores zero, and so does a point whose mean distances to the rest of
    its own group and to the nearest other group are both zero.
    """
    _, codes, counts = np.unique(labels, return_inverse=True, return_counts=True)
    points = np.arange(len(codes))
    own = codes == np.arange(len(counts))[:, None]
    # The sum of the distances from each point (column) to each group (row).
    sums = np.stack([distances[:, members].sum(axis=1) for members in own])
    # A point is not one of the others in its own group.
    means = sums / np.maximum(counts[:, None] - own, 1)
    within = means[codes, points]
    nearest = np.where(own, np.inf, means).min(axis=0)
    spread = np.maximum(within, nearest)
    return np.divide(
        nearest - within,
        spread,
        out=np.zeros(len(codes)),
        where=(counts[codes] > 1) & (spread > 0),
    )


def log_silhouettes(basename, ordered, unclustered, names, distances):
    """Log the mean silhouette of the clustered sequences, and of each cluster.

    ordered is as ordered_clusters returns it, and names and distances as
    patristic_distances returns them.  The Unique sequences are not a cluster,
    so they are left out of the score.
    """
    if len(ordered) < 2:
        logger.info(
            "No silhouette: it needs at least two clusters (found %d)",
            len(ordered),
        )
        return
    column = {name: index for index, name in enumerate(names)}
    members = [column[seqid] for seqids, _clade in ordered for seqid in seqids]
    labels = np.repeat(np.arange(len(ordered)), [len(s) for s, _clade in ordered])
    scores = silhouette(distances[np.ix_(members, members)], labels)
    logger.info(
        "Mean silhouette %.3f over the %d sequences of %d clusters, "
        "leaving out %d Unique sequences",
        scores.mean(),
        len(members),
        len(ordered),
        len(unclustered),
    )
    for index in range(len(ordered)):
        logger.info(
            "Silhouette of %s.%d: %.3f", basename, index, scores[labels == index].mean()
        )


def write_clusters(basename, seq_idx, tree, ordered, unclustered):
    """Write output files: clusters & unique as FASTA, tree as phyloXML.

    The files are named after basename, the base name of the input file.
    seq_idx maps each identifier to its record, as read_sequences returns it,
    and ordered lists the clusters as ordered_clusters returns them, in the
    order of their numbers.
    """

    def write_cluster(seqids, fname):
        """Write the sequences of the given identifiers to a FASTA file."""
        records = [seq_idx[seqid] for seqid in seqids]
        write_fasta(records, fname)
        logger.info("Wrote %s (%d sequences)", fname, len(records))

    for i, (seqids, clade) in enumerate(ordered):
        write_cluster(seqids, f"{basename}.{i}")
        # Evenly spaced hues, one per cluster.
        rgb = colorsys.hsv_to_rgb(i / len(ordered), 0.7, 0.9)
        clade.color = BranchColor(*(int(x * 255) for x in rgb))
        clade.width = 2
    if unclustered:
        write_cluster(sorted(unclustered), f"{basename}.Unique")

    treefname = f"{basename}.xml"
    Phylo.write(tree, treefname, "phyloxml")
    logger.info("Wrote %s", treefname)
