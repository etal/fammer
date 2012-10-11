#!/usr/bin/env python

"""Use TMalign for multiple structure alignment.

Procedure:

1. Run all-vs-all pairwise alignments w/ TMalign; keep the alignments and also
   the TM-score (a metric of structural similarity/alignability).
2. Build a graph of structure similarities based on TM-scores. Determine the
   minimum spanning tree (MST) for the graph.
3. Take the pairwise alignments corresponding to the edges (connected nodes) of
   the MST, and use them as seeds in MAFFT (--localpair).
4. Remove duplicated sequences and all-gap columns from the resulting alignment.


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


def align_structs(pdb_fnames, seed_fnames=None):
    """Align multiple PDB structures using TM-align.

    Returns a list of aligned SeqRecords.
    """
    # 1. Align all-v-all structure pairs with TM-align
    allpairs = []
    for idx, ref_pdbfn in enumerate(pdb_fnames):
        assert ' ' not in ref_pdbfn
        for eqv_pdbfn in pdb_fnames[idx+1:]:
            assert eqv_pdbfn != ref_pdbfn
            logging.info("Aligning %s to %s", eqv_pdbfn, ref_pdbfn)
            tm_output = subprocess.check_output(['TMalign', ref_pdbfn, eqv_pdbfn])
            tm_seqpair = read_tmalign_as_seqrec_pair(tm_output, ref_pdbfn, eqv_pdbfn)
            allpairs.append(tm_seqpair)

    # In case of 2 structs, no need to combine alignments -- we're done
    if len(allpairs) == 1:
        recs = allpairs[0][:2]
        return alnutils.remove_empty_cols(recs)

    # 2. Resolve MST pairs & write seed tempfiles
    tmp_seed_fnames = []
    for seedpair in mst_pairs(allpairs):
        # fd, seedfn = tempfile.mkstemp(text=True)
        # SeqIO.write(seedpair, seedfn, 'fasta')
        # SeqIO.write(seedpair, os.fdopen(fd), 'fasta')
        with tempfile.NamedTemporaryFile('w+', delete=False) as handle:
            SeqIO.write(seedpair, handle, 'fasta')
            tmp_seed_fnames.append(handle.name)

    # 3. Use MAFFT to combine TMalign'd pairs into a multiple alignment;
    seq_fd, seq_fname = tempfile.mkstemp(text=True)
    # Create a blank file to appease MAFFT
    os.write(seq_fd, '')
    os.close(seq_fd)
    mafft_output = subprocess.check_output(['mafft',
        '--quiet', '--amino', '--localpair',
        '--maxiterate', '1000']
        + list(itertools.chain(*[('--seed', sfn)
                                 for sfn in (seed_fnames or []) + tmp_seed_fnames]))
        + [seq_fname])
    # Clean up
    os.remove(seq_fname)
    for sfn in tmp_seed_fnames:
        os.remove(sfn)

    # 4. Emit the aligned sequences
    recs = SeqIO.parse(StringIO(mafft_output), 'fasta')
    recs = clean_and_dedupe_seqs(recs)
    return list(alnutils.remove_empty_cols(recs))


def read_tmalign_as_seqrec_pair(tm_output, ref_id, eqv_id):
    """Create a pair of SeqRecords from TMalign output."""
    lines = tm_output.splitlines()
    # Extract the TM-score (measure of structure similarity)
    # Take the mean of the (two) given TM-scores -- not sure which is reference
    tmscores = []
    for line in lines:
        if line.startswith('TM-score'):
            tmscores.append(float(line.split(None, 2)[1]))
    tmscore = sum(tmscores) / len(tmscores)
    # Extract the sequence alignment
    lastlines = lines[-5:]
    assert lastlines[0].startswith('(":"') # (":" denotes the residues pairs
    assert not lastlines[-1].strip()
    refseq, eqvseq = lastlines[1].strip(), lastlines[3].strip()
    return (SeqRecord(Seq(refseq), id=ref_id,
                      description="TMalign TM-score=%f" % tmscore),
            SeqRecord(Seq(eqvseq), id=eqv_id,
                      description="TMalign TM-score=%f" % tmscore),
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


def tmscore_from_description(text):
    for token in text.split():
        if token.startswith('TM-score'):
            return float(token.split('=', 1)[1])


def clean_and_dedupe_seqs(records, best_score=False):
    """Remove the _seed_ prefix and omit duplicated records."""
    if best_score:
        seen = {}
    else:
        seen = set()
    for record in records:
        # Remove the _seed_ prefix from each sequence ID
        if record.id.startswith('_seed_'):
            record.id = record.id[len('_seed_'):]
            record.name = record.id
        # Check for duplicates.
        if best_score:
            # If a previously seen PDB was aligned better (per the TM-score),
            # defer to that one
            tmscore = tmscore_from_description(record.description)
            if record.id in seen and seen[record.id] >= tmscore:
                # This PDB was aligned better previously; skip
                continue
            seen[record.id] = tmscore
        else:
            # Keep a duplicate sequence if it was aligned differently
            ident = (record.id, str(record.seq))
            if ident in seen:
                continue
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

