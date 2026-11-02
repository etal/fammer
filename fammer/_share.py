"Helpers shared across multiple Fammer sub-commands."

import collections
import tempfile
from os.path import getsize, join

from Bio import SeqIO

from .tasks import ext, sh


class UsageError(ValueError):
    """An error in a command's input, reported without a traceback.

    It is a ValueError, so that a library caller catching ValueError for a bad
    argument still catches it.
    """


# HMMer utilities

# The e-value at or below which a profile's HMM counts as hitting a sequence,
# when a set of sequences is scored against one profile.
MEMBER_EVALUE = 0.001


def get_hmm_hits(hmmfname, seqfname, evalue):
    """Search a sequence database with an HMM profile database.

    Return a mapping of hit sequence names to the best profile, with the
    e-value and bit score of its best domain, as parse_scanned_hits does.

    The sequence file is read as FASTA. A file of zero bytes has no hits, and
    hmmsearch is not run on it, since it exits with an error on one; a file of
    only whitespace has no hits either, and a misformatted file still fails.
    """
    if getsize(seqfname) == 0:
        return {}
    with tempfile.NamedTemporaryFile(mode="w+", suffix=".tbl") as tbl:
        sh(
            "hmmsearch --tformat fasta --noali --notextw "
            f"-E {evalue} --domE {evalue} --tblout {tbl.name} "
            f"{hmmfname} {seqfname} > /dev/null"
        )
        tbl.seek(0)
        return parse_scanned_hits(tbl)


def get_hmm_evalues(hmmfname, seqfname):
    """Map each sequence that an HMM profile hits to its best-domain e-value."""
    hits = get_hmm_hits(hmmfname, seqfname, MEMBER_EVALUE)
    return {seqname: evalue for seqname, (_prof, evalue, _score) in hits.items()}


def get_hmm_scores(hmmfname, seqfname):
    """Map each sequence that an HMM profile hits to its best-domain bit score."""
    hits = get_hmm_hits(hmmfname, seqfname, MEMBER_EVALUE)
    return {seqname: score for seqname, (_prof, _evalue, score) in hits.items()}


def build_hmm(records, hmmfname):
    """Build an HMM profile from aligned records, as fammer build does.

    The records are written as Stockholm beside the HMM file, then hmmbuild
    reads them with its default options.
    """
    stkfname = ext(hmmfname, "stk")
    SeqIO.write(records, stkfname, "stockholm")
    sh(f"hmmbuild {hmmfname} {stkfname}")


def leave_one_out_hmms(rows, held_out, workdir):
    """Build, for each row to hold out, an HMM from the other aligned rows.

    The rows are the aligned records of one profile, and held_out the
    identifiers of the rows to hold out in turn; the other rows stay in every
    HMM. Generates tuples of (identifier, HMM file name) in row order. Each HMM
    is written to the same file in workdir, so a caller uses it before taking
    the next one.
    """
    hmmfname = join(workdir, "held-out.hmm")
    for index, row in enumerate(rows):
        if row.id not in held_out:
            continue
        build_hmm(rows[:index] + rows[index + 1 :], hmmfname)
        yield row.id, hmmfname


def parse_scanned_hits(tblfile):
    """Parse a file stream from hmmsearch --tblout.

    Returns a dict of {target_seq_name: (best_profile_name, e-value, bit score)}

    Each sequence is matched to the profile whose best single domain has the
    smallest e-value, the first such row on a tie. The e-value and the bit
    score are those of that domain, from the same row.

    Example contents::

        #                                                                  --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
        # target name           accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
        #   ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
        CDPK-consensus          -          Testfam              -           1.9e-153  508.4  25.4  2.1e-153  508.2  17.6   1.0   1   0   0   1   1   1   1 -
        gi|254575025|pdb|3HX4|A -          Testfam              -           1.4e-150  498.9   5.7  1.6e-150  498.7   4.0   1.0   1   0   0   1   1   1   1 {<Apicomplexa(e)>}Chain A, Crystal Structure Of Cdpk1 Of Toxoplasma Gondii, Tgme49_1014 Presence Of Calcium
        TAIR|locus:2011201      -          Testfam              -           2.4e-137  455.2   0.0  2.7e-137  455.1   0.0   1.0   1   0   0   1   1   1   1 symbol:CDPK1 species:3702 "Arabidopsis thaliana"  AT1G18890
    """  # noqa: E501
    best = {}
    for line in tblfile:
        if line.startswith("#"):
            continue
        (
            seqname,
            _acc1,
            profname,
            _acc2,
            seq_eval,
            seq_score,
            seq_bias,
            dom_eval,
            dom_score,
            dom_bias,
            _etc,
        ) = line.split(None, 10)
        if seqname in best and float(dom_eval) >= best[seqname][1]:
            continue
        best[seqname] = (profname, float(dom_eval), float(dom_score))
    return best


# Helpers


def repeated_names(names):
    """The names that occur more than once, sorted."""
    counts = collections.Counter(names)
    return sorted(name for name, count in counts.items() if count > 1)


def ungapped(seq):
    """The sequence as a string without gap characters, in upper case."""
    return str(seq).replace("-", "").replace(".", "").upper()


def write_fasta(records, dest):
    """Write records to a path or handle as FASTA without wrapping lines."""
    SeqIO.write(records, dest, "fasta-2line")
