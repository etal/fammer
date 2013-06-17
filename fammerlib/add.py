"Add sequences from a scanned database to the profile tree."

import logging
import os
import tempfile
from os.path import basename, isdir, isfile, join

from Bio import SeqIO

from .tasks import ext, noext, sh
from ._share import get_hmm_evalues, parse_scanned_hits


def cmd_add(args):
    # Validation
    assert isdir(args.basedir), "Not a directory: %s" % args.basedir
    hmmdb = args.basedir.rstrip('/') + '_all.hmm'
    assert isfile(hmmdb), "Didn't find " + hmmdb
    assert isfile(args.sequence_db)
    assert args.evalue
    assert 0 < args.min_length <= 1, \
            "Minimum hit length (-L, --min-length) must be between 0 and 1."

    # Traverse the directory tree -- if there's an unbuilt HMM that we'll need,
    # best to identify it now and raise an error
    fastas_and_hmms = list(walk_profiles(args.basedir))
    # Scan & classify the target database, just like in cmd_scan
    dbhits = get_hmm_data(hmmdb, args.sequence_db, args.evalue)
    assert dbhits
    index = SeqIO.index(args.sequence_db, 'fasta')
    assert index

    for fasta, hmm in fastas_and_hmms:
        logging.info("Processing %s ...", noext(fasta))
        # HMMscan the corresponding FASTA file; record the max (worst) e-value
        seq_evals = get_hmm_evalues(hmm, fasta)
        max_eval = max(seq_evals.values())

        # Select the hits matching the current profile
        hits = list(filter_profile(dbhits, basename(noext(fasta))))
        # if not hits:
            # logging.info("No hits for profile %s", hmm)
            # continue
        # Filter hits by this profile's maximum (worst) e-value
        hits = list(filter_evalue(hits, max_eval))
        if not hits:
            # logging.info("No hits for profile %s with e-value <= %s",
            #              hmm, max_eval)
            continue

        # Align the hit sequences (retrieved from the target DB)
        aln = get_hmm_alignment(hmm, hits, index)

        # Calculate minimum sequence length based on HMM profile length
        hmm_len = get_hmm_len(hmm)
        min_len_aa = int(hmm_len * args.min_length)

        # Fetch the aligned sequence regions as SeqRecords
        hit_seqs = []
        for record in aln:
            # Read the domain hit start/end position from Stockholm record ID
            acc, _offset = record.id.split('/')
            start, end = map(int, _offset.split('-'))
            # Filter hits by minimum length
            hit_len = end - start + 1   # HMMer uses 1-based numbering
            if hit_len < min_len_aa:
                continue
            # Remove the "/123-456" suffix from aligned record keys
            record.id = acc
            # Remove gaps from the hit sequence
            record.seq._data = str(record.seq).\
                    replace('-', '').replace('.', '').upper()
            hit_seqs.append(record)

        # ENH: extend hits that are close to the edge of the profile (e.g. 5aa) to
        # fill sequence to the edges (no gaps)
        #   e.g. if hmmstart = 4, add another 3aa sliced from the full sequence
        #   to the aligned bit -- if only 2aa available, then take that much
        # - can I tell from the alignment? (# terminal gaps)

        # Add the accepted hits to the FASTA profile
        logging.info("Adding %d hits to profile %s", len(hit_seqs), fasta)
        with open(fasta, 'a') as fafile:
            SeqIO.write(hit_seqs, fafile, 'fasta')


def walk_profiles(topdir):
    """Iterate through paired FASTA and HMM profiles in the directory tree."""
    for dirpath, _subdirs, fnames in os.walk(topdir):
        fasta_fnames = [join(dirpath, fn)
                        for fn in fnames if fn.endswith('.fasta')]
        for fasta in fasta_fnames:
            hmm_fname = ext(fasta, 'hmm')
            if not isfile(hmm_fname):
                raise ValueError("Unbuild .hmm for %s" % fasta)
            yield fasta, hmm_fname


def get_hmm_data(hmmfname, seqfname, evalue):
    """Search a sequence database with an HMM profile.

    Return a mapping of hit sequence names to the best profile and e-value.
    """
    with tempfile.NamedTemporaryFile(suffix='.tbl') as tbl:
        sh('hmmsearch --noali --notextw '
           '-E %s --domE %s --tblout %s %s %s > /dev/null'
           % (evalue, evalue, tbl.name, hmmfname, seqfname))
        tbl.seek(0)
        # Parse the table; get the domain score for each sequence
        hits = parse_scanned_hits(tbl)
    return hits


def get_hmm_alignment(hmmfname, hits, index):
    fullseqs = [index[seqname] for seqname, _e in hits]
    assert fullseqs
    tmp_seqfname = hmmfname + '.fa'
    with open(tmp_seqfname, 'w+') as handle:
        SeqIO.write(fullseqs, handle, 'fasta')
        handle.flush()
    with tempfile.NamedTemporaryFile(suffix='.stk') as alnfile:
        sh(('hmmsearch --noali -A %s %s %s > /dev/null')
            % (alnfile.name, hmmfname, tmp_seqfname))
        alnfile.seek(0)
        alignment = list(SeqIO.parse(alnfile, 'stockholm'))
    os.remove(tmp_seqfname)
    return alignment


def get_hmm_len(hmmfname):
    """Read the HMM profile length from an HMM file."""
    with open(hmmfname) as hmmfile:
        for line in hmmfile:
            if line.startswith('LENG'):
                key, length = line.split()
                assert key == 'LENG'
                return int(length)


def filter_profile(hits, profile_name):
    """Select the hits matching the given profile name.

    Hits are a dict of: {sequence name: (best profile, e-value)}

    Generates tuples of: (sequence name, e-value)
    """
    for seqname, data in hits.iteritems():
        profname, evalue = data
        if profname == profile_name:
            yield (seqname, evalue)


def filter_evalue(hits, max_evalue):
    """Skip hits that are above the given e-value."""
    for hit in hits:
        if hit[1] <= max_evalue:
            yield hit



