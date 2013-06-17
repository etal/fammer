"Functions shared across multiple Fammer sub-commands."

import collections
import tempfile

from Bio.File import as_handle

from .tasks import sh


# HMMer utilities

def get_hmm_evalues(hmmfname, seqfname):
    with tempfile.NamedTemporaryFile(suffix='.tbl') as tbl:
        sh(('hmmsearch --noali --notextw -E 0.001 --domE 0.001 '
            '--tblout %s %s %s > /dev/null')
            % (tbl.name, hmmfname, seqfname))
        tbl.seek(0)
        # Parse the table; get the domain score for each sequence
        hits = parse_scanned_hits(tbl)
    return dict((seqname, prof_eval[1])
                for seqname, prof_eval in hits.iteritems())


def parse_scanned_hits(tblfile):
    """Parse a file stream from hmmsearch --tblout.

    Returns a dict of {target_seq_name: (best_profile_name, e-value)}

    Example contents:

#                                                                  --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name           accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#   ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
CDPK-consensus          -          Testfam              -           1.9e-153  508.4  25.4  2.1e-153  508.2  17.6   1.0   1   0   0   1   1   1   1 -
gi|254575025|pdb|3HX4|A -          Testfam              -           1.4e-150  498.9   5.7  1.6e-150  498.7   4.0   1.0   1   0   0   1   1   1   1 {<Apicomplexa(e)>}Chain A, Crystal Structure Of Cdpk1 Of Toxoplasma Gondii, Tgme49_1014 Presence Of Calcium
TAIR|locus:2011201      -          Testfam              -           2.4e-137  455.2   0.0  2.7e-137  455.1   0.0   1.0   1   0   0   1   1   1   1 symbol:CDPK1 species:3702 "Arabidopsis thaliana"  AT1G18890

    """
    topevals = collections.defaultdict(int)
    topprofs = {}
    _lines = 0
    for line in tblfile:
        if line.startswith('#'):
            continue
        _lines += 1
        (seqname, _acc1, profname, _acc2,
                seq_eval, seq_score, seq_bias,
                dom_eval, dom_score, dom_bias,
                _etc) = line.split(None, 10)
        # Match each sequence to its highest-scoring profile
        if seqname in topevals and float(dom_eval) >= topevals[seqname]:
            continue
        topevals[seqname] = float(dom_eval)
        topprofs[seqname] = profname
    return dict((seq, (prof, topevals[seq]))
                 for seq, prof in topprofs.iteritems())


# Helpers

# def touch(task):
#     """Just create a file if it's not already present."""
#     with open(task.target, 'w'):
#         pass


def write_fasta(records, fname):
    """Write a FASTA file without wrapping lines."""
    with as_handle(fname, 'w+') as outfile:
        for rec in records:
            descr = rec.description.strip()
            if descr:
                outfile.write(">%s %s\n%s\n" % (rec.id, descr, rec.seq))
            else:
                outfile.write(">%s\n%s\n" % (rec.id, rec.seq))

