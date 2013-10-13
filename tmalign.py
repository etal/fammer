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

"""

import argparse
import sys
import logging

from Bio import SeqIO

from fammerlib.tmalign import align_structs


if __name__ == '__main__':
    AP = argparse.ArgumentParser(__doc__)
    AP.add_argument('pdbfnames', nargs='+',
            help="PDB file names of the structures to align.")
    AP.add_argument('-s', '--seed', action='append', default=[],
            help="Aligned FASTA sequences to include in the output alignment.")
    AP.add_argument('-o', '--output',
                  default=sys.stdout,
                  type=argparse.FileType('w+'),
                  help="Filename for the output FASTA alignment.")
    args = AP.parse_args()
    logging.basicConfig(level=logging.INFO,
            format="%(module)s [@%(lineno)s]: %(message)s")

    if args.seed:
        print >>sys.stderr, "tmalign: Using seed alignments", \
                ' '.join(args.seed)
        out_seqrecs = align_structs(args.pdbfnames, seed_fnames=args.seed)
    else:
        out_seqrecs = align_structs(args.pdbfnames)
    SeqIO.write(out_seqrecs, args.output, 'fasta')
    print >>sys.stderr, "tmalign: Wrote", args.output.name
