#!/usr/bin/env python2.7

"""A toolkit for protein superfamily sequence profile alignment and curation.

Manages a hierarchical directory tree containing sequence sets (FASTA files)
corresponding to each subfamily of a protein superfamily.
"""
from __future__ import absolute_import

import logging
import argparse

from fammerlib.build import cmd_build
from fammerlib.scan import cmd_scan
from fammerlib.add import cmd_add
from fammerlib.refine import cmd_refine
from fammerlib.cluster import cmd_cluster
from fammerlib.update_fasta import cmd_update_fasta


if __name__ == '__main__':
    AP = argparse.ArgumentParser(
            description=__doc__,
            epilog="Contact Eric <etal@uga.edu> for help.")

    # Global options
    AP.add_argument('-q', '--quiet',
            action='store_true',
            help="Don't print status messages, only warnings and errors.")
    AP_subparsers = AP.add_subparsers(
            help='sub-commands (use with -h for more info)')

    # Sub-command: build
    P_build = AP_subparsers.add_parser('build',
            help='Build profiles from a given directory tree.')
    P_build.add_argument('basedir',
            help='Top of the directory tree where the sequences sets are.')
    P_build.add_argument('--align-script', dest='align_script', metavar='FILE',
            help='Script to align a FASTA file, writing Clustal format.')
    P_build.add_argument('--clean',
            action='store_true',
            help='Clean auxilliary files after building.')
    P_build.add_argument('--no-pdb',
            dest='pdb', action='store_false', default=True,
            help="Don't use PDB structures when constructing alignments.")
    P_build.add_argument('--hmmer',
            action='store_true',
            help='Build profiles for HMMer.')
    P_build.add_argument('--mapgaps',
            action='store_true',
            help='Build profiles for MAPGAPS.')
    P_build.add_argument('--tree',
            action='store_true',
            help='Write a Newick tree indicating the directory structure.')
    P_build.set_defaults(func=cmd_build)

    # Sub-command: scan
    P_scan = AP_subparsers.add_parser('scan',
            help='Scan and classify input sequences using a set of profiles.')
    P_scan.add_argument('profile_db',
            help='HMM profile name (or basename).')
    P_scan.add_argument('target',
            help='Sequences to scan, in unaligned FASTA format.')
    P_scan.add_argument('-E', '--evalue',
            default='0.0001',
            help='E-value cutoff for significant search hits.')
    P_scan.add_argument('--table',
            action='store_true',
            help='Print a table of sequence-to-profile assignments.')
    P_scan.add_argument('-o', '--seqout',
            help='Write a FASTA file of all sequences matching any profile.')
    P_scan.add_argument('-O', '--seqsets',
            action='store_true',
            help='Write FASTA files of the sequences matching each profile.')
    P_scan.add_argument('-A', '--align',
            action='store_true',
            help='Make a FASTA alignment of sequences matching each profile.')
    P_scan.add_argument('-n', '--include',
            help='Only match profiles whose name contains this string(s).')
    P_scan.add_argument('-x', '--exclude',
            help='Exclude matches to profiles whose name contains this string.')
    P_scan.set_defaults(func=cmd_scan)

    # Sub-command: add
    P_add = AP_subparsers.add_parser('add',
            help="""Scan a target database with the given HMM profile set. Add
            hits that meet acceptance thresholds to the profile FASTA files.""")
    P_add.add_argument('basedir',
            help="""Top of the directory tree where the sequences sets are.
            There must also be an HMM profile database at <basedir>_all.hmm.""")
    P_add.add_argument('sequence_db',
            help='Sequence database to scan for new additions to the profiles.')
    P_add.add_argument('-E', '--evalue',
            default='1e-20',
            help='E-value cutoff for significant search hits.')
    P_add.add_argument('-L', '--min-length',
            default=0.8, type=float,
            help="""Minimum length, as a proportion of the profile length, of
            hits to accept into the profiles.""")
    # ENH: -u, --unclassified : exclude? put parent-level hits here?
    P_add.set_defaults(func=cmd_add)

    # Sub-command: refine
    P_refine = AP_subparsers.add_parser('refine',
            help='Leave-one-out validation of HMM profiles.')
    P_refine.add_argument('basedir',
            help='Top of the directory tree where the sequences sets are.')
    P_refine.add_argument('-n', '--dry-run',
            action='store_true',
            help='Do a "dry run"; don\'t modify the original sequence sets.')
    P_refine.add_argument('-u', '--unclassified',
            help="""The set of unclassified/unassigned sequences (FASTA file).
            By default, this is <basedir>/<basedir>-Unique.fasta.""")
    P_refine.add_argument('-x', '--exclude',
            help="Exclude this profile from evaluation.")
    P_refine.set_defaults(func=cmd_refine)

    # Sub-command: update-fasta
    P_update = AP_subparsers.add_parser('update-fasta', 
            help="""Replace original FASTA sequence sets with the ungapped
            sequences from the corresponding alignment (.aln) files, sorted by
            decreasing length.""")
    P_update.add_argument('basedir',
            help='Top of the directory tree where the sequences sets are.')
    P_update.set_defaults(func=cmd_update_fasta)

    # Sub-command: cluster
    P_cluster = AP_subparsers.add_parser('cluster',
            help="""Split a sequence set into smaller clusters.""")
    P_cluster.add_argument('sequences',
            help='Clustal alignment (.aln) or unaligned FASTA (.fasta).')
    P_cluster.set_defaults(func=cmd_cluster)

    args = AP.parse_args()
    # Handle global options here
    if args.quiet:
        logging.basicConfig(level=logging.WARNING,
                format="%(module)s: %(message)s")
    else:
        logging.basicConfig(level=logging.INFO,
                format="%(module)s [@%(lineno)s]: %(message)s")
    # Pass the rest along to the sub-command implementation
    args.func(args)

