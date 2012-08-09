#!/usr/bin/env python

"""Use TMalign for multiple structure alignment.

Run all-vs-all pairwise alignments w/ TMalign; keep the alignments and also the
metric of similarity/distance.

Build a graph of similarities between structs. Determine the minimum spanning
tree (MST) for the graph. Take the pairwise alignments corresponding to the
edges (connected nodes) of the MST, use them as seeds in MAFFT (--localpair).

In degenerate cases (1 or 2 structs), skip the graph part.

Depends:
    - TMalign
    - MAFFT
    - networkx
    - Biopython
    - biofrills
"""

import itertools
import logging
import os
import subprocess
import tempfile
from cStringIO import StringIO

import networkx

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biofrills import alnutils


def align_structs(pdb_fnames):
    # Align each other PDB against the "reference" PDB
    allpairs = []
    # 1. Align all-v-all pairs
    for idx, ref_pdbfn in enumerate(pdb_fnames):
        assert ' ' not in ref_pdbfn
        for eqv_pdbfn in pdb_fnames[idx+1:]:
            assert eqv_pdbfn != ref_pdbfn
            logging.info("Aligning %s to %s", eqv_pdbfn, ref_pdbfn)
            tm_output = subprocess.check_output(['TMalign', ref_pdbfn, eqv_pdbfn])
            # Convert the TMalign output into a FASTA sequence pair
            tm_seqpair = read_tmalign_as_seqrec_pair(tm_output, ref_pdbfn, eqv_pdbfn)
            allpairs.append(tm_seqpair)
            # pairfafname = eqv_pdbfn + '.tm.fa'
            # SeqIO.write(tm_seqpair, pairfafname, 'fasta')
            # seedfnames.append(pairfafname)

    # 2. Resolve MST pairs & write seed tempfiles
    seedfnames = []
    for seedpair in mst_pairs(allpairs):
        # seedpair is a tuple of 2 SeqRecords
        seedfn = tempfile.mkstemp(text=True)[1]
        SeqIO.write(list(seedpair), seedfn, 'fasta')
        seedfnames.append(seedfn)

    # 3. Use MAFFT to combine TMalign'd pairs into a multiple alignment;
    seq_fname = tempfile.mkstemp(text=True)[1]
    with open(seq_fname, 'w') as handle:
        # Create a blank file to appease MAFFT
        pass
    mafft_output = subprocess.check_output(['mafft',
        '--quiet', '--amino', '--localpair',
        '--maxiterate', '1000']
        + list(itertools.chain(*[('--seed', sfn) for sfn in seedfnames]))
        + [seq_fname])
    # Clean up
    os.remove(seq_fname)
    for sfn in seedfnames:
        os.remove(sfn)

    # Emit the aligned sequences
    recs = SeqIO.parse(StringIO(mafft_output), 'fasta')
    recs = clean_and_dedupe_seqs(recs)
    return alnutils.remove_empty_cols(recs)


def read_tmalign_as_seqrec_pair(tm_output, ref_id, eqv_id):
    """Create a pair of SeqRecords from TMalign output."""
    lines = tm_output.splitlines()
    # TODO - extract TM-score; larger means a better fit
    # Extract the TM-score (measure of structure similarity)
    for line in lines:
        if line.startswith('Aligned'):
            for token in line.split():
                if token.startswith('TM-score='):
                    tmscore = float(token[:-1].split('=')[1])
                    break
            break
    # Extract the sequence alignment
    lastlines = lines[-5:]
    assert lastlines[0].startswith('(":"') # (":" denotes the residues pairs
    assert not lastlines[-1].strip()
    refseq, eqvseq = lastlines[1].strip(), lastlines[3].strip()
    return (SeqRecord(Seq(refseq), id=ref_id, description="TMalign"),
            SeqRecord(Seq(eqvseq), id=eqv_id, description="TMalign"),
            tmscore)


def mst_pairs(pairs):
    """Given all pairwise distances, determine the minimal spanning subset.

    Convert pairwise distances to an undirected graph, determine the
    minumum spanning tree, and emit the minimal list of edges to connect all
    nodes.

    Input: iterable of (SeqRecord, SeqRecord, distance)
    Output: iterable of (SeqRecord, SeqRecord)
    """
    G = networkx.Graph()
    for left, right, score in pairs:
        G.add_edge(left, right, weight=1.0/score)
    mst = networkx.minimum_spanning_edges(G, data=False)
    return list(mst)


def clean_and_dedupe_seqs(records):
    """Remove the _seed_ prefix and omit duplicated records."""
    seen = set()
    for record in records:
        # Remove the _seed_ prefix from each sequence ID
        if record.id.startswith('_seed_'):
            record.id = record.id[len('_seed_'):]
        # Skip exact duplicates. The same PDB can be aligned differently to
        # other PDBs, so also check the aligned sequence.
        ident = (record.id, str(record.seq))
        if ident not in seen:
            seen.add(ident)
            yield record


# --- command line ------------------------------------------------------

def process_args(args):
    out_seqrecs = align_structs(args.pdbfnames)
    SeqIO.write(out_seqrecs, args.output, 'fasta')


if __name__ == '__main__':
    import argparse
    import sys
    try:
        from esbglib.sugar import log_config
        log_config()
    except ImportError:
        pass

    AP = argparse.ArgumentParser(__doc__)
    AP.add_argument('pdbfnames', nargs='+',
            help="PDB file names of the structures to align.")
    AP.add_argument('-o', '--output',
                  default=sys.stdout,
                  type=argparse.FileType('w+'),
                  help="Filename for the output FASTA alignment.")
    process_args(AP.parse_args())

