"Add sequences from a scanned database to the profile tree."

import logging
import os
import tempfile
from contextlib import closing
from os.path import basename, isdir, isfile, join

from Bio import SeqIO
from Bio.Seq import Seq

from ._share import (
    UsageError,
    get_hmm_evalues,
    get_hmm_hits,
    ungapped,
    write_fasta,
)
from .tasks import ext, noext, sh

logger = logging.getLogger(__name__)


def cmd_add(args):
    if not isdir(args.basedir):
        raise UsageError(f"Not a directory: {args.basedir}")
    hmmdb = args.basedir.rstrip("/") + "_all.hmm"
    if not isfile(hmmdb):
        raise UsageError(f"No HMM profile database at {hmmdb}")
    if not isfile(args.sequence_db):
        raise UsageError(f"No sequence file at {args.sequence_db}")
    if not 0 < args.min_length <= 1:
        raise UsageError(
            "Minimum hit length (-L, --min-length) must be between 0 and 1."
        )

    # Traverse the directory tree -- if there's an unbuilt HMM that we'll need,
    # best to identify it now and raise an error
    fastas_and_hmms = list(walk_profiles(args.basedir))
    # Scan & classify the target database, just like in cmd_scan
    dbhits = get_hmm_hits(hmmdb, args.sequence_db, args.evalue)
    if not dbhits:
        raise UsageError(f"No profile hits in {args.sequence_db}")
    with closing(SeqIO.index(args.sequence_db, "fasta")) as index:
        for fasta, hmm in fastas_and_hmms:
            add_hits(fasta, hmm, dbhits, index, args.min_length)


def add_hits(fasta, hmm, dbhits, index, min_length):
    """Append the database hits that a profile accepts to its FASTA file."""
    logger.info("Processing %s ...", noext(fasta))
    # Select the hits matching the current profile
    hits = list(filter_profile(dbhits, basename(noext(fasta))))
    if not hits:
        return
    # Skip the sequences the profile already holds, so that a repeated run
    # adds each sequence once
    hits = list(filter_present(hits, member_names(fasta)))
    if not hits:
        return

    # Score the profile's own sequences with its own HMM; the worst e-value is
    # the acceptance threshold.
    seq_evals = get_hmm_evalues(hmm, fasta)
    if not seq_evals:
        logger.warning(
            "Skipping %s: its HMM profile hits none of its own sequences at "
            "e-value 0.001 or better, so there is no threshold to accept a hit",
            fasta,
        )
        return
    max_eval = max(seq_evals.values())
    hits = list(filter_evalue(hits, max_eval))
    if not hits:
        return

    # Align the hit sequences (retrieved from the target DB)
    aln = get_hmm_alignment(hmm, hits, index)

    # Calculate minimum sequence length based on HMM profile length
    hmm_len = get_hmm_len(hmm)
    min_len_aa = int(hmm_len * min_length)

    hit_seqs = []
    for acc, record in longest_domains(aln, min_len_aa):
        # Remove the "/123-456" suffix from aligned record keys
        record.id = acc
        # Remove gaps from the hit sequence. The per-residue posterior
        # probabilities from hmmsearch -A do not fit the ungapped
        # sequence, so drop them first.
        record.letter_annotations = {}
        record.seq = Seq(ungapped(record.seq))
        hit_seqs.append(record)
    if not hit_seqs:
        return

    # ENH: extend hits that are close to the edge of the profile (e.g. 5aa) to
    # fill sequence to the edges (no gaps)
    #   e.g. if hmmstart = 4, add another 3aa sliced from the full sequence
    #   to the aligned bit -- if only 2aa available, then take that much
    # - can I tell from the alignment? (# terminal gaps)

    # Add the accepted hits to the FASTA profile
    logger.info("Adding %d hits to profile %s", len(hit_seqs), fasta)
    with open(fasta, "a") as fafile:
        write_fasta(hit_seqs, fafile)


def walk_profiles(topdir):
    """Iterate through paired FASTA and HMM profiles in the directory tree."""
    for dirpath, _subdirs, fnames in os.walk(topdir):
        fasta_fnames = [join(dirpath, fn) for fn in fnames if fn.endswith(".fasta")]
        for fasta in fasta_fnames:
            hmm_fname = ext(fasta, "hmm")
            if not isfile(hmm_fname):
                raise UsageError(
                    f"No .hmm profile for {fasta}; run fammer build --hmmer first"
                )
            yield fasta, hmm_fname


def get_hmm_alignment(hmmfname, hits, index):
    """Align the hit sequences to the HMM profile with hmmsearch -A.

    Returns the Stockholm records, one per aligned domain, each named
    ``<sequence name>/<start>-<end>``.
    """
    fullseqs = [index[seqname] for seqname, _e in hits]
    with (
        tempfile.NamedTemporaryFile(mode="w+", suffix=".fa") as seqfile,
        tempfile.NamedTemporaryFile(mode="w+", suffix=".stk") as alnfile,
    ):
        write_fasta(fullseqs, seqfile)
        seqfile.flush()
        sh(f"hmmsearch --noali -A {alnfile.name} {hmmfname} {seqfile.name} > /dev/null")
        alnfile.seek(0)
        return list(SeqIO.parse(alnfile, "stockholm"))


def get_hmm_len(hmmfname):
    """Read the HMM profile length from an HMM file."""
    with open(hmmfname) as hmmfile:
        for line in hmmfile:
            if line.startswith("LENG"):
                key, length = line.split()
                assert key == "LENG"
                return int(length)


def filter_profile(hits, profile_name):
    """Select the hits matching the given profile name.

    Hits are a dict of: {sequence name: (best profile, e-value)}

    Generates tuples of: (sequence name, e-value)
    """
    for seqname, data in hits.items():
        profname, evalue = data
        if profname == profile_name:
            yield (seqname, evalue)


def member_names(fasta):
    """The identifiers of the sequences in a profile's FASTA file."""
    with SeqIO.parse(fasta, "fasta") as records:
        return {record.id for record in records}


def filter_present(hits, members):
    """Skip the hits whose sequence name is among the given members."""
    for hit in hits:
        if hit[0] not in members:
            yield hit


def longest_domains(records, min_len_aa):
    """Select the longest domain of each sequence that hmmsearch -A aligned.

    The records are named ``<sequence name>/<start>-<end>``, one per aligned
    domain, so a sequence with several domains appears several times. Domains
    shorter than min_len_aa are skipped; of two equally long domains, the
    first is kept.

    Generates tuples of: (sequence name, record)
    """
    longest = {}
    for record in records:
        # Sequence names may themselves contain "/"
        acc, offset = record.id.rsplit("/", 1)
        start, end = map(int, offset.split("-"))
        hit_len = end - start + 1  # HMMer uses 1-based numbering
        if hit_len < min_len_aa:
            continue
        if acc in longest and longest[acc][0] >= hit_len:
            continue
        longest[acc] = (hit_len, record)
    for acc, (_hit_len, record) in longest.items():
        yield acc, record


def filter_evalue(hits, max_evalue):
    """Skip hits that are above the given e-value."""
    for hit in hits:
        if hit[1] <= max_evalue:
            yield hit
