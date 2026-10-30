"Helpers shared across multiple Fammer sub-commands."

import collections
import tempfile

from Bio import SeqIO

from .tasks import sh


class UsageError(ValueError):
    """An error in a command's input, reported without a traceback.

    It is a ValueError, so that a library caller catching ValueError for a bad
    argument still catches it.
    """


# HMMer utilities


def get_hmm_hits(hmmfname, seqfname, evalue):
    """Search a sequence database with an HMM profile database.

    Return a mapping of hit sequence names to the best profile and e-value.
    """
    with tempfile.NamedTemporaryFile(mode="w+", suffix=".tbl") as tbl:
        sh(
            "hmmsearch --noali --notextw "
            f"-E {evalue} --domE {evalue} --tblout {tbl.name} "
            f"{hmmfname} {seqfname} > /dev/null"
        )
        tbl.seek(0)
        return parse_scanned_hits(tbl)


def get_hmm_evalues(hmmfname, seqfname):
    """Map each sequence that an HMM profile hits to its best e-value."""
    hits = get_hmm_hits(hmmfname, seqfname, 0.001)
    return {seqname: evalue for seqname, (_prof, evalue) in hits.items()}


def parse_scanned_hits(tblfile):
    """Parse a file stream from hmmsearch --tblout.

    Returns a dict of {target_seq_name: (best_profile_name, e-value)}

    Example contents::

        #                                                                  --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
        # target name           accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
        #   ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
        CDPK-consensus          -          Testfam              -           1.9e-153  508.4  25.4  2.1e-153  508.2  17.6   1.0   1   0   0   1   1   1   1 -
        gi|254575025|pdb|3HX4|A -          Testfam              -           1.4e-150  498.9   5.7  1.6e-150  498.7   4.0   1.0   1   0   0   1   1   1   1 {<Apicomplexa(e)>}Chain A, Crystal Structure Of Cdpk1 Of Toxoplasma Gondii, Tgme49_1014 Presence Of Calcium
        TAIR|locus:2011201      -          Testfam              -           2.4e-137  455.2   0.0  2.7e-137  455.1   0.0   1.0   1   0   0   1   1   1   1 symbol:CDPK1 species:3702 "Arabidopsis thaliana"  AT1G18890
    """  # noqa: E501
    topevals = collections.defaultdict(int)
    topprofs = {}
    _lines = 0
    for line in tblfile:
        if line.startswith("#"):
            continue
        _lines += 1
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
        # Match each sequence to its highest-scoring profile
        if seqname in topevals and float(dom_eval) >= topevals[seqname]:
            continue
        topevals[seqname] = float(dom_eval)
        topprofs[seqname] = profname
    return dict((seq, (prof, topevals[seq])) for seq, prof in topprofs.items())


# Helpers


def ungapped(seq):
    """The sequence as a string without gap characters, in upper case."""
    return str(seq).replace("-", "").replace(".", "").upper()


def write_fasta(records, dest):
    """Write records to a path or handle as FASTA without wrapping lines."""
    SeqIO.write(records, dest, "fasta-2line")
