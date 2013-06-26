"""Extract phylogenetically-based clusters from a sequence set.

Uses FastTree to quickly build a tree with branch support values, then extracts
well-supported clades from the tree and writes the corresponding sequence sets
to FASTA files. Unclustered sequences are written to another "Unique" file.
Also writes a phyloXML tree file (.xml) showing clusters as colorized clades.

Depends:
    biopython (pip install biopython)
    fasttree (apt-get install fasttree)
    mafft (apt-get install mafft)

"""
from __future__ import division

import logging
import os.path
import subprocess
import tempfile
from cStringIO import StringIO

from Bio import AlignIO, SeqIO, Phylo
from Bio.Phylo.BaseTree import BranchColor
from Bio.Graphics.ColorSpiral import ColorSpiral

from biofrills import alnutils


def cmd_cluster(args):
    """Cut on longest branches, then group by remaining supported clusters."""
    # Args:
    #     sequences: Target sequence file (foo.aln, clustal format)
    # Opts (ENH):
    #     --size: MIN_CLUSTER_SIZE (e.g. 4)
    #     --confidence: MIN_CONFIDENCE (e.g. 0.7, 0.9) -- needed?
    #     --tree: precomputed tree, incl. supports -- good idea?
    #     --outgroup: needed?
    assert os.path.isfile(args.sequences)
    tree = load_tree(args.sequences)

    # Extract clusters of acceptable size
    MIN_CLUSTER_SIZE = max(4, int(round(tree.count_terminals() ** .4)))
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
    for clade in tree.find_clades(order='level'):
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

    write_clusters(args.sequences, tree, clusters, unclustered)


def load_tree(seqfname):
    """Load an alignment, build & prep a tree, return the tree object."""
    if seqfname.endswith('.aln'):
        aln = AlignIO.read(seqfname, 'clustal')
    elif seqfname.endswith('.fasta'):
        # Run MAFFT quickly
        alndata = subprocess.check_output(['mafft', '--quiet', '--auto',
                                           seqfname])
        aln = AlignIO.read(StringIO(alndata), 'fasta')
    else:
        raise ValueError("Input sequences must be a Clustal alignment (.aln) "
                         "or unaligned FASTA (.fasta)")

    # Use conserved (less-gappy) blocks to build the tree
    aln = alnutils.blocks(aln, 0.4)
    with tempfile.NamedTemporaryFile(mode='w') as tmp:
        AlignIO.write(aln, tmp, 'fasta')
        tmp.flush()
        treedata = subprocess.check_output(['fasttree',
                                            '-pseudo', '-gamma', '-wag',
                                            tmp.name])
    tree = Phylo.read(StringIO(treedata), 'newick')

    # Collapse weakly supported splits
    confs = [c.confidence
             for c in tree.find_clades()
             if c.confidence is not None]
    # ENH: accept min_confidence as an option
    min_confidence = sum(confs) / len(confs)
    tree.collapse_all(lambda c: c.confidence < min_confidence)
    tree.ladderize(reverse=True)
    tree.root.branch_length = 0.0
    return tree


def get_clades_to_cut(tree):
    """Select long internal branches to cut the tree at."""
    # Index internal branches by length & clade
    tip_branch_lengths = []
    inner_branch_lengths = []
    for clade in tree.find_clades():
        if clade.is_terminal():
            tip_branch_lengths.append(clade.branch_length)
        else:
            inner_branch_lengths.append(clade.branch_length)
    # tip_mean = sum(tip_branch_lengths) / len(tip_branch_lengths)
    tip_mean = geom_mean(tip_branch_lengths)
    clades_to_cut = [cl for cl in tree.get_nonterminals()
                    if cl.branch_length >= tip_mean]
    clades_to_cut.sort(key=lambda cl: cl.branch_length, reverse=True)
    logging.info("# to cut (up to): %s", len(clades_to_cut))
    return clades_to_cut


def geom_mean(arr):
    """Geometric mean of an array of numbers."""
    prod = reduce(float.__mul__, map(float, arr))
    return pow(prod, 1./len(arr))

   
def write_clusters(seqfname, tree, clusters, unclustered):
    """Write output files: clusters & unique as FASTA, tree as phyloXML."""
    is_aln = seqfname.endswith('.aln')
    seq_idx = SeqIO.to_dict(SeqIO.parse(seqfname,
                                        'clustal' if is_aln else 'fasta'))
    def write_cluster(cluster, fname):
        """Write the sequences of cluster tips to a FASTA file."""
        records = [seq_idx[seqid] for seqid in sorted(cluster)]
        with open(fname, 'w+') as handle:
            for rec in records:
                write_fasta(rec, handle, do_ungap=is_aln)
        logging.info("Wrote %s (%d sequences)", fname, len(records))

    colors = [BranchColor(*map(lambda x: int(x*255), rgb))
            for rgb in ColorSpiral().get_colors(len(clusters))]
    for i, item in enumerate(sorted(clusters.iteritems(), reverse=True,
                                    key=lambda kv: len(kv[1]))):
        clade, cluster = item
        write_cluster(cluster, os.path.basename(seqfname) + '.' + str(i))
        clade.color = colors[i]
        clade.width = 2
    if unclustered:
        write_cluster(unclustered, os.path.basename(seqfname) + '.Unique')

    treefname = os.path.basename(seqfname) + '.xml'
    Phylo.write(tree, treefname, 'phyloxml') 
    logging.info("Wrote %s", treefname)


def write_fasta(record, handle, do_ungap=False):
    """Ungap and write a SeqRecord to a handle in single-line FASTA format."""
    if record.description:
        topline = '>%s %s' % (record.id, record.description)
    else:
        topline = '>' + record.id + '\n'
    seq = str(record.seq).strip()
    if do_ungap:
        seq = seq.replace('-', '').replace('.', '').upper()
    handle.write(topline + '\n' + seq + '\n')

