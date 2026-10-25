"""Command-line entry points: the fammer and fammer-tmalign programs."""

import logging
import argparse
import sys

from Bio import SeqIO

from . import tmalign
from .build import cmd_build
from .scan import cmd_scan
from .add import cmd_add
from .refine import cmd_refine
from .cluster import cmd_cluster
from .update_fasta import cmd_update_fasta


def _configure_logging(quiet=False):
    """Set the log format and level shared by both programs."""
    if quiet:
        logging.basicConfig(level=logging.WARNING, format="%(module)s: %(message)s")
    else:
        logging.basicConfig(
            level=logging.INFO, format="%(module)s [@%(lineno)s]: %(message)s"
        )


def main():
    """Run the fammer program."""
    AP = argparse.ArgumentParser(
        description="""A toolkit for protein superfamily sequence profile
            alignment and curation. Manages a hierarchical directory tree
            containing sequence sets (FASTA files) corresponding to each
            subfamily of a protein superfamily.""",
        epilog="Contact Eric <etal@uga.edu> for help.",
    )

    # Global options
    AP.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Don't print status messages, only warnings and errors.",
    )
    AP_subparsers = AP.add_subparsers(
        required=True, help="sub-commands (use with -h for more info)"
    )

    # Sub-command: build
    P_build = AP_subparsers.add_parser(
        "build", help="Build profiles from a given directory tree."
    )
    P_build.add_argument(
        "basedir", help="Top of the directory tree where the sequences sets are."
    )
    P_build.add_argument(
        "--align-script",
        dest="align_script",
        metavar="FILE",
        help="Script to align a FASTA file, writing Clustal format.",
    )
    P_build.add_argument(
        "--clean", action="store_true", help="Clean auxilliary files after building."
    )
    P_build.add_argument(
        "--no-pdb",
        dest="pdb",
        action="store_false",
        default=True,
        help="Don't use PDB structures when constructing alignments.",
    )
    P_build.add_argument(
        "--hmmer", action="store_true", help="Build profiles for HMMer."
    )
    P_build.add_argument(
        "--mapgaps", action="store_true", help="Build profiles for MAPGAPS."
    )
    P_build.add_argument(
        "--tree",
        action="store_true",
        help="Write a Newick tree indicating the directory structure.",
    )
    P_build.set_defaults(func=cmd_build)

    # Sub-command: scan
    P_scan = AP_subparsers.add_parser(
        "scan", help="Scan and classify input sequences using a set of profiles."
    )
    P_scan.add_argument("profile_db", help="HMM profile name (or basename).")
    P_scan.add_argument("target", help="Sequences to scan, in unaligned FASTA format.")
    P_scan.add_argument(
        "-E",
        "--evalue",
        default="0.0001",
        help="E-value cutoff for significant search hits.",
    )
    P_scan.add_argument(
        "--table",
        action="store_true",
        help="Print a table of sequence-to-profile assignments.",
    )
    P_scan.add_argument(
        "-o",
        "--seqout",
        help="Write a FASTA file of all sequences matching any profile.",
    )
    P_scan.add_argument(
        "-O",
        "--seqsets",
        action="store_true",
        help="Write FASTA files of the sequences matching each profile.",
    )
    P_scan.add_argument(
        "-A",
        "--align",
        action="store_true",
        help="Make a FASTA alignment of sequences matching each profile.",
    )
    P_scan.add_argument(
        "-n",
        "--include",
        help="Only match profiles whose name contains this string(s).",
    )
    P_scan.add_argument(
        "-x",
        "--exclude",
        help="Exclude matches to profiles whose name contains this string.",
    )
    P_scan.set_defaults(func=cmd_scan)

    # Sub-command: add
    P_add = AP_subparsers.add_parser(
        "add",
        help="""Scan a target database with the given HMM profile set. Add
            hits that meet acceptance thresholds to the profile FASTA files.""",
    )
    P_add.add_argument(
        "basedir",
        help="""Top of the directory tree where the sequences sets are.
            There must also be an HMM profile database at <basedir>_all.hmm.""",
    )
    P_add.add_argument(
        "sequence_db",
        help="Sequence database to scan for new additions to the profiles.",
    )
    P_add.add_argument(
        "-E",
        "--evalue",
        default="1e-20",
        help="E-value cutoff for significant search hits.",
    )
    P_add.add_argument(
        "-L",
        "--min-length",
        default=0.8,
        type=float,
        help="""Minimum length, as a proportion of the profile length, of
            hits to accept into the profiles.""",
    )
    # ENH: -u, --unclassified : exclude? put parent-level hits here?
    P_add.set_defaults(func=cmd_add)

    # Sub-command: refine
    P_refine = AP_subparsers.add_parser(
        "refine", help="Leave-one-out validation of HMM profiles."
    )
    P_refine.add_argument(
        "basedir", help="Top of the directory tree where the sequences sets are."
    )
    P_refine.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help='Do a "dry run"; don\'t modify the original sequence sets.',
    )
    P_refine.add_argument(
        "-u",
        "--unclassified",
        help="""The set of unclassified/unassigned sequences (FASTA file).
            By default, this is <basedir>/<basedir>-Unique.fasta.""",
    )
    P_refine.add_argument(
        "-x", "--exclude", help="Exclude this profile from evaluation."
    )
    P_refine.set_defaults(func=cmd_refine)

    # Sub-command: update-fasta
    P_update = AP_subparsers.add_parser(
        "update-fasta",
        help="""Replace original FASTA sequence sets with the ungapped
            sequences from the corresponding alignment (.aln) files, sorted by
            decreasing length.""",
    )
    P_update.add_argument(
        "basedir", help="Top of the directory tree where the sequences sets are."
    )
    P_update.set_defaults(func=cmd_update_fasta)

    # Sub-command: cluster
    P_cluster = AP_subparsers.add_parser(
        "cluster", help="""Split a sequence set into smaller clusters."""
    )
    P_cluster.add_argument(
        "sequences", help="Clustal alignment (.aln) or unaligned FASTA (.fasta)."
    )
    P_cluster.set_defaults(func=cmd_cluster)

    args = AP.parse_args()
    # Handle global options here
    _configure_logging(args.quiet)
    # Pass the rest along to the sub-command implementation
    args.func(args)


def tmalign_main():
    """Run the fammer-tmalign program."""
    AP = argparse.ArgumentParser(
        description=tmalign.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    AP.add_argument(
        "pdbfnames", nargs="+", help="PDB file names of the structures to align."
    )
    AP.add_argument(
        "-s",
        "--seed",
        action="append",
        default=[],
        help="Aligned FASTA sequences to include in the output alignment.",
    )
    AP.add_argument(
        "-o",
        "--output",
        default=sys.stdout,
        type=argparse.FileType("w+"),
        help="Filename for the output FASTA alignment.",
    )
    args = AP.parse_args()
    _configure_logging()

    if args.seed:
        print("tmalign: Using seed alignments", " ".join(args.seed), file=sys.stderr)
    out_seqrecs = tmalign.align_structs(args.pdbfnames, seed_fnames=args.seed)
    SeqIO.write(out_seqrecs, args.output, "fasta")
    print("tmalign: Wrote", args.output.name, file=sys.stderr)
