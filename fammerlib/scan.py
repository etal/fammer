"Scan a sequence set using a HMMer 3.0 profile database."
# ENH - use the built HMM/MAPGAP & Newick tree to scan&classify seqs

import collections
import logging
from os.path import basename, isdir, isfile

from Bio import SeqIO

from .tasks import sh
from ._share import parse_scanned_hits


def cmd_scan(args):
    """Use HMMer to classify sequences like MAPGAPS does.

    Target sequence identifiers must be unique!
    """
    # Verify: HMM profile exists
    if isdir(args.profile_db):
        hmmdb = args.profile_db.rstrip('/') + '_all.hmm'
    else:
        hmmdb = args.profile_db
    assert isfile(hmmdb), "Didn't find " + hmmdb
    assert isfile(args.target)

    # ENH?: remove table ([target].tbl) with option --clean
    tblfname = basename(args.target) + '.tbl'
    sh(('hmmsearch -E %s --domE %s --noali ' \
        '--tblout %s %s %s > /dev/null')
        % (args.evalue, args.evalue, tblfname, hmmdb, args.target))
    # Get the best hit for each sequence
    with open(tblfname) as tblfile:
        scanhits = parse_scanned_hits(tblfile)
    if not scanhits:
        logging.warn("No hits in %s", args.target)
        return

    # Filter hits for include/exclude strings
    if args.include:
        includes = args.include.split(',')
        for key, val in scanhits.items():
            profname, _evalue = val
            if not any(incl in profname for incl in includes):
                del scanhits[key]
    if args.exclude:
        excludes = args.exclude.split(',')
        for key, val in scanhits.items():
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
    profcounts = collections.Counter(
        [prof for prof, score in scanhits.itervalues()])
    max_width = max(map(len, profcounts.iterkeys()))
    for profname, count in sorted(profcounts.iteritems(),
            # sort by decreasing count, then alphabetically
            key=lambda pair: (-pair[1], pair[0])):
        print profname.ljust(max_width), ':', count


def write_table(scanhits):
    """Table of accessions & best-matching profile.

    Rows of:
    [sequence id] [tab] [profile name]
    """
    # ENH: include e-value cutoff, sequence description
    # [sequence id] \t [profile name] \t [e-value] \t [description...]
    for seqname, data in sorted(scanhits.iteritems()):
        profname, score = data
        print "%s\t%s" % (seqname, profname)


def write_seqout(scanhits, sourcefname, outfname):
    """Write a single FASTA file containing all hits to any profiles."""
    seq_idx = SeqIO.index(sourcefname, 'fasta')
    with open(outfname, 'w+') as outfile:
        for acc in scanhits:
            block = seq_idx.get_raw(acc)
            outfile.write(block)
    logging.info("Wrote %s", outfname)


def write_seqsets(scanhits, sourcefname, align=False):
    """Write FASTA files of the seqs matching each profile."""
    seq_idx = SeqIO.index(sourcefname, 'fasta')
    # NB: scanhits looks like {seq_id: (profile_name, score)}
    # Flip to: {profile_name: [sequence ids...]}
    profiles = collections.defaultdict(list)
    for seqname, prof_score in scanhits.iteritems():
        profiles[prof_score[0]].append(seqname)
    # For each entry there, fetch all seq ids from the file & write
    for profname, seqnames in sorted(profiles.iteritems()):
        outfname = "%s.%s.fasta" % (basename(sourcefname), profname)
        with open(outfname, 'w+') as outfile:
            for rec in seqnames:
                block = seq_idx.get_raw(rec)
                outfile.write(block)
        logging.info("Wrote %s", outfname)

    # 3. hmmalign the written sequence sets (if requested)
    # if align: ...
    # ENH - use the built HMM/MAPGAP & Newick tree to scan&classify&align seqs
    # HMMer:
    # do cmd_scan -O to get each seq's family
    #   ungap along the way
    # with NamedTemporaryFile(suffix='.hmm'):
    #     hmmfetch -o tmp.hmm [hmm_db] [family_name]
    #     hmmalign tmp.hmm [tgt_fam] > tgt_result.sto
    # SeqIO.convert('tgt_result.sto', 'stockholm', 'tgt_result.seq', 'fasta)

