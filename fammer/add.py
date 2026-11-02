"Add sequences from a scanned database to the profile tree."

import functools
import logging
import os
import tempfile
from collections import defaultdict
from contextlib import closing
from os.path import basename, isdir, isfile, join

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ._share import (
    MEMBER_EVALUE,
    UsageError,
    best_hits,
    get_hmm_rows,
    get_hmm_scores,
    is_catch_all,
    leave_one_out_hmms,
    ungapped,
    write_fasta,
)
from .tasks import ends_with_newline, ext, noext, sh

logger = logging.getLogger(__name__)

# A profile of at most this many sequences takes its acceptance threshold from
# leave-one-out scores; see profile_threshold.
LEAVE_ONE_OUT_MAX = 10


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
    profiles = {}
    for fasta, hmm in walk_profiles(args.basedir):
        name = basename(noext(fasta))
        if name in profiles:
            raise UsageError(
                f"Two profiles are named {name}: {profiles[name][0]} and {fasta}"
            )
        profiles[name] = (fasta, hmm)
    # A sequence already in any profile is not added again, to that profile or
    # to another, so a repeated run adds each sequence once
    present = set().union(*(member_names(fasta) for fasta, _hmm in profiles.values()))
    # Scan the target database with every profile of the tree, as cmd_scan
    # does, but assign each sequence only among the profiles that have a FASTA
    # file to receive it
    rows = get_hmm_rows(hmmdb, args.sequence_db, args.evalue)
    if not rows:
        raise UsageError(f"No profile hits in {args.sequence_db}")
    rows = [row for row in filter_present(rows, present) if row[1] in profiles]
    threshold = functools.cache(lambda name: profile_threshold(*profiles[name]))
    dbhits = assign_hits(rows, threshold)
    with closing(SeqIO.index(args.sequence_db, "fasta")) as index:
        for name, (fasta, hmm) in profiles.items():
            logger.info("Processing %s ...", noext(fasta))
            hits = list(filter_profile(dbhits, name))
            if hits:
                add_hits(fasta, hmm, hits, threshold(name), index, args.min_length)


def assign_hits(rows, threshold):
    """Assign each sequence of hmmsearch table rows to one profile.

    The rows are tuples of (sequence name, profile name, e-value, bit score),
    and threshold maps a profile name to its acceptance threshold in bits, or
    to None if it has none. A sequence goes to the profile whose best domain
    has the smallest e-value, as best_hits selects it, unless that profile is
    a catch-all: then it goes to the subfamily profile with the smallest
    e-value among those whose threshold its bit score reaches, if there is
    one, since a catch-all holds only the sequences that fit no subfamily.

    Returns a dict of {sequence name: (profile name, e-value, bit score)}.
    """
    rows_of = defaultdict(list)
    for row in rows:
        rows_of[row[0]].append(row)
    best = best_hits(rows)
    for seqname, (profname, _evalue, _score) in best.items():
        if not is_catch_all(profname):
            continue
        accepted = [
            row
            for row in rows_of[seqname]
            if not is_catch_all(row[1])
            and threshold(row[1]) is not None
            and row[3] >= threshold(row[1])
        ]
        if accepted:
            best[seqname] = best_hits(accepted)[seqname]
    return best


def add_hits(fasta, hmm, hits, min_score, index, min_length):
    """Append the hits that a profile accepts to its FASTA file.

    The hits are tuples of (sequence name, bit score), and min_score is the
    profile's acceptance threshold (see profile_threshold), or None.
    """
    # The weakest of the profile's own sequences sets the acceptance
    # threshold, in bits. A bit score, unlike an e-value, does not depend on
    # the number of sequences searched, so the profile's own scores and the
    # target's are comparable.
    if min_score is None:
        logger.warning(
            "Skipping %s: its HMM profile hits none of its own sequences at "
            "e-value %s or better, so there is no threshold to accept a hit",
            fasta,
            MEMBER_EVALUE,
        )
        return
    hits = list(filter_score(hits, min_score))
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

    # Add the accepted hits to the FASTA profile. A profile whose last line
    # has no line ending gets one first, so the first header starts a line.
    logger.info("Adding %d hits to profile %s", len(hit_seqs), fasta)
    needs_newline = not ends_with_newline(fasta)
    with open(fasta, "a") as fafile:
        if needs_newline:
            fafile.write("\n")
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


def profile_threshold(fasta, hmm):
    """The lowest best-domain bit score among a profile's own sequences.

    The sequences are the rows of the profile's alignment (.aln), from which
    fammer build made its HMM, that are also in its FASTA file, scored without
    gaps; a profile without an alignment uses its FASTA file. A record that
    add appended since the last build is not a row, so it cannot lower the
    threshold before the HMM is rebuilt with it.

    A profile of at most LEAVE_ONE_OUT_MAX sequences scores each against an HMM
    built from the alignment without its row, since against an HMM built
    partly from itself a sequence of a small profile scores far higher than an
    unseen member of the family does. A larger profile scores them against its
    own HMM. Only the sequences hit at MEMBER_EVALUE or better count; without
    any, the result is None.
    """
    alnfname = ext(fasta, "aln")
    if not isfile(alnfname):
        return min(get_hmm_scores(hmm, fasta).values(), default=None)
    names = member_names(fasta)
    with SeqIO.parse(alnfname, "clustal") as parsed:
        rows = list(parsed)
    members = [
        SeqRecord(Seq(ungapped(row.seq)), id=row.id, description="")
        for row in rows
        if row.id in names
    ]
    with tempfile.TemporaryDirectory(prefix="fammer-") as workdir:
        seqfname = join(workdir, "members.fasta")
        write_fasta(members, seqfname)
        if len(rows) < 2 or len(members) > LEAVE_ONE_OUT_MAX:
            scores = get_hmm_scores(hmm, seqfname)
        else:
            scores = {}
            held_out = {member.id for member in members}
            for name, lean in leave_one_out_hmms(rows, held_out, workdir):
                score = get_hmm_scores(lean, seqfname).get(name)
                if score is not None:
                    scores[name] = score
    return min(scores.values(), default=None)


def get_hmm_alignment(hmmfname, hits, index):
    """Align the hit sequences to the HMM profile with hmmsearch -A.

    Returns the Stockholm records, one per aligned domain, each named
    ``<sequence name>/<start>-<end>``.
    """
    fullseqs = [index[hit[0]] for hit in hits]
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

    Hits are a dict of: {sequence name: (best profile, e-value, bit score)}

    Generates tuples of: (sequence name, bit score)
    """
    for seqname, (profname, _evalue, score) in hits.items():
        if profname == profile_name:
            yield (seqname, score)


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


def filter_score(hits, min_score):
    """Skip hits whose bit score is below the given one."""
    for hit in hits:
        if hit[1] >= min_score:
            yield hit
