"Scan a sequence set using a HMMer 3.0 profile database."
# ENH - use the built HMM/MAPGAP & Newick tree to scan&classify seqs

import collections
import logging
from contextlib import closing
from os.path import basename, isdir, isfile

from Bio import SeqIO

from ._share import UsageError, get_hmm_hits

logger = logging.getLogger(__name__)


def cmd_scan(args):
    """Use HMMer to classify sequences like MAPGAPS does.

    Target sequence identifiers must be unique!
    """
    # Verify: HMM profile exists
    if isdir(args.profile_db):
        hmmdb = args.profile_db.rstrip("/") + "_all.hmm"
    else:
        hmmdb = args.profile_db
    if not isfile(hmmdb):
        raise UsageError(f"No HMM profile database at {hmmdb}")
    if not isfile(args.target):
        raise UsageError(f"No sequence file at {args.target}")

    # Get the best hit for each sequence
    scanhits = get_hmm_hits(hmmdb, args.target, args.evalue)
    if not scanhits:
        logger.warning("No hits in %s", args.target)
        return

    # Filter hits for include/exclude strings
    if args.include:
        includes = args.include.split(",")
        for key, val in list(scanhits.items()):
            profname, _evalue = val
            if not any(incl in profname for incl in includes):
                del scanhits[key]
    if args.exclude:
        excludes = args.exclude.split(",")
        for key, val in list(scanhits.items()):
            profname, _evalue = val
            if any(excl in profname for excl in excludes):
                del scanhits[key]

    # Write a table of just the best hit for each sequence
    if args.table:
        write_table(scanhits)
    else:
        write_summary(scanhits)
    if args.seqout:
        write_seqout(scanhits, args.target, args.seqout)
    if args.seqsets:
        write_seqsets(scanhits, args.target)


def write_summary(scanhits):
    """Summary of # hits per profile.

    Rows of:
    [profile name] : [#hits]
    """
    profcounts = collections.Counter([prof for prof, _score in scanhits.values()])
    max_width = max(map(len, profcounts.keys()))
    for profname, count in sorted(
        profcounts.items(),
        # sort by decreasing count, then alphabetically
        key=lambda pair: (-pair[1], pair[0]),
    ):
        print(profname.ljust(max_width), ":", count)


def write_table(scanhits):
    """Table of accessions & best-matching profile.

    Rows of:
    [sequence id] [tab] [profile name]
    """
    # ENH: include e-value cutoff, sequence description
    # [sequence id] \t [profile name] \t [e-value] \t [description...]
    for seqname, data in sorted(scanhits.items()):
        profname, _score = data
        print(f"{seqname}\t{profname}")


def write_seqout(scanhits, sourcefname, outfname):
    """Write a single FASTA file containing all hits to any profiles."""
    with (
        closing(SeqIO.index(sourcefname, "fasta")) as seq_idx,
        open(outfname, "wb+") as outfile,
    ):
        for acc in scanhits:
            outfile.write(seq_idx.get_raw(acc))
    logger.info("Wrote %s", outfname)


def write_seqsets(scanhits, sourcefname):
    """Write FASTA files of the seqs matching each profile."""
    # NB: scanhits looks like {seq_id: (profile_name, score)}
    # Flip to: {profile_name: [sequence ids...]}
    profiles = collections.defaultdict(list)
    for seqname, prof_score in scanhits.items():
        profiles[prof_score[0]].append(seqname)
    # For each entry there, fetch all seq ids from the file & write
    with closing(SeqIO.index(sourcefname, "fasta")) as seq_idx:
        for profname, seqnames in sorted(profiles.items()):
            outfname = f"{basename(sourcefname)}.{profname}.fasta"
            with open(outfname, "wb+") as outfile:
                for rec in seqnames:
                    outfile.write(seq_idx.get_raw(rec))
            logger.info("Wrote %s", outfname)
