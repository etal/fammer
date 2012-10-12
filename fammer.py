#!/usr/bin/env python2.7

"""A toolkit for protein superfamily sequence profile alignment and curation.

Manages a hierarchical directory tree containing sequence sets (FASTA files)
corresponding to each subfamily of a protein superfamily.
"""

import collections
import itertools
import logging
import os
import subprocess
import sys
import tempfile
import time
from glob import glob
from os.path import abspath, basename, isdir, isfile, join

from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biocma import biocma
from biofrills import consensus, alnutils

from tasks import Task, ext, noext, sh, is_empty
import tmalign


# === Build ===========================================

class Result(object):
    """Output of taskify_subdirs. Collection of tasks.

    Contains references to the tasks for the main alignment, the PDB seed
    alignment, and, optionally, the corresponding HMMer .hmm profile and/or
    MAPGAPS .cma and .tpl.
    """
    def __init__(self, aln, pdbseq=None, hmm=None, cma=None, tpl=None):
        self.aln = aln
        self.pdbseq = pdbseq
        self.hmm = hmm
        self.cma = cma
        self.tpl = tpl

    def __str__(self):
        return str(self.aln)
    # ENH: __cmp__ or __eq__ and __le__ and full_ordering
    # ENH: build() and clean() methods?


def cmd_build(args):
    # Verify: 'base' directory exists
    assert isdir(args.basedir)
    base = args.basedir.rstrip('/')

    # Write the directory layout as a Newick tree
    if args.tree:
        tree_str = treeify_subdirs(base) + ';\n'
        with open(base + '.nwk', 'w+') as outfile:
            outfile.write(tree_str)
        logging.info("Wrote %s.nwk", base)

    # Build all alignments, plus HMM .hmm or MAPGAPS .cma/.tpl
    base_result = taskify_subdirs(base, args.hmmer, args.mapgaps, args.pdb, 0)

    # HMM profiles
    if args.hmmer:
        T_base_all_hmm = Task(base + '_all.hmm',
                action=all_hmm,
                depends=base_result.hmm,
                )
        T_base_all_hmm.build()
        # Clean everything
        if args.clean:
            T_base_all_hmm.clean()

    # MAPGAPS profiles
    if args.mapgaps:
        T_mg_run_map = Task(base + '.mpa',
                action=mg_run_map,
                depends=[base_result.cma, base_result.tpl])
        T_mg_run_map.build()
        if args.clean:
            T_mg_run_map.clean()

    if not args.hmmer and not args.mapgaps:
        # Just build the alignments
        base_result.aln.build()
        if args.clean:
            base_result.aln.clean()


# To build the HMM (and CMA) profiles

def taskify_subdirs(topdir, hmmer, mapgaps, use_pdb, level):
    """Recursively define tasks for a directory tree.

    Return a task to build the HMM profile for topdir.

    this.hmm
        this.aln
            this.fasta -- consenses from each group/family
                this/subdir1.hmm (recurse)
                this/subdir2.hmm (recurse)
                this/...
                this/[others.hmm]
                --> task: hmmemit from each of deez

    How do we know which sub-hmms to build?
        subdirs of this/ --> group names --> [this/groupname.hmm]
        this/*.fasta --> [this/family.hmm]

    If mapgaps is True, build CMA profiles and templates too.
    """
    this = topdir.rstrip('/')
    subdirs = filter(isdir, 
            [join(topdir, sd) for sd in sorted(os.listdir(topdir))])
    # Each subtask pair: (HMM, CMA)-making tasks for the subdir
    # Groups with their own families within -- recurse
    subtask_group_results = [taskify_subdirs(sd, hmmer, mapgaps, use_pdb, level+1)
                             for sd in sorted(subdirs)]
    # Families / tips of the profile tree -- build from scratch
    subtask_family_results = []
    subfamily_fastas = glob(join(topdir, '*.fasta'))
    for subfa in sorted(subfamily_fastas):
        # Skip the group sequence sets, they're already covered
        if noext(subfa) in subdirs:
            continue
        subresult = Result(Task(ext(subfa, 'aln'),
                action=align_fasta,
                depends=subfa,
                cleans=ext(subfa, 'seq')))
        if hmmer:
            subresult.hmm = Task(ext(subfa, 'hmm'),
                    action=aln2hmm,
                    depends=subresult.aln,
                    cleans=ext(subfa, 'stk'))
        if mapgaps:
            subresult.cma = Task(ext(subfa, 'cma'),
                    action=mg_aln2cma,
                    kwargs={'level': level+1},
                    depends=subresult.aln,
                    cleans=[ext(subfa, 'cons.cma'),
                            ext(subfa, 'cons_iron.cma')])
        subtask_family_results.append(subresult)

    subtask_family_results.sort(key=str)    # Needed? We sort FASTAs above

    # Structural alignment of PDBs in this dir; reuse subfamily PDB alignments
    these_pdbs = glob(join(topdir, '*.pdb'))
    sub_pdb_seqs = [sgr.pdbseq for sgr in subtask_group_results]
    task_pdbseq = Task(this + '.pdb.seq',
            action=align_pdbs,
            kwargs={'use_pdb': use_pdb},
            depends=these_pdbs + sub_pdb_seqs)

    # Aggregate those profile consensus sequences & make a meta-profile
    result = Result(Task(this + '.aln',
                         action=align_profiles,
                         kwargs={'use_pdb': use_pdb},
                         depends=[r.aln
                                  for r in (subtask_group_results +
                                            subtask_family_results)
                                 ] + [task_pdbseq],
                         cleans=ext(map(str, subtask_group_results +
                                        subtask_family_results),
                                    'cons.seq') + [this + '.families.fa',
                                                   this + '.families.seq',
                                                   this + '.seq']),
                    pdbseq=task_pdbseq)
    if hmmer:
        result.hmm = Task(this + '.hmm',
                action=aln2hmm,
                depends=([result.aln] +
                    [sgr.hmm for sgr in subtask_group_results] +
                    [sfr.hmm for sfr in subtask_family_results]),
                cleans=this + '.stk')
    if mapgaps:
        result.tpl = Task(this + '.tpl',
                action=mg_aln2cma,
                kwargs={'level': level},
                depends=result.aln,
                cleans=[this + '.cons.cma', this + '.cons_iron.cma'])
        result.cma = Task(this + '.cma',
                action=mg_cat_cma,
                depends=list(itertools.chain(*[
                    (sgr.tpl, sgr.cma)
                    for sgr in subtask_group_results]
                    )) + [sfr.cma for sfr in subtask_family_results])

    return result


# Actions

# def touch(task):
#     """Just create a file if it's not already present."""
#     with open(task.target, 'w'):
#         pass


def align_fasta(task):
    """Align a FASTA file with MAFFT. Clustal output."""
    seq = ext(task.depends[0], 'seq')
    sh("mafft --quiet --amino --reorder --maxiterate 1000 "
       "--genafpair --ep 0.123 %s > %s"
       % (task.depends[0], seq))
    # XXX fast version
    # sh("mafft --quiet --amino --reorder --auto %s > %s" % (task.depends[0], seq))
    # Convert FASTA to "pressed" (single-row) Clustal
    records = list(SeqIO.parse(seq, 'fasta'))
    # check for 'X' characters in the sequences -- these cause problems
    for rec in records:
        if 'X' in str(rec.seq):
            logging.warn('Sequence %r contains unknown residue X' % rec.id)
    max_id_len = max(len(r.id) for r in records)
    with open(task.target, 'w+') as outfile:
        outfile.write('CLUSTAL X (-like) multiple sequence alignment\n\n')
        outfile.writelines(
                ['%s %s\n' % (rec.id.ljust(max_id_len), rec.seq)
                 for rec in records])


def align_pdbs(task, sub_pdb_seqs=(), use_pdb=None):
    """Create a structure-based sequence alignment from PDB files."""
    if not use_pdb:
        # Just touch the '.pdb.seq' file; don't use TM-align
        with open(task.target, 'a'):
            return

    pdbs = []
    sub_pdb_seqs = []
    for elem in task.depends:
        if str(elem).endswith('.pdb'):
            pdbs.append(elem)
        else:
            sub_pdb_seqs.append(str(elem))
    # Scan existing PDB alignments to choose a reference PDB from each
    # sub_pdb_seqs = filter(isfile, sub_pdb_seqs)
    for sub_pdb_fname in map(str, sub_pdb_seqs):
        if is_empty(sub_pdb_fname):
            # Dummy seed alignment -- look in that dir for a .pdb
            #   ENH - recursively
            sub_pdbs = glob(join(sub_pdb_fname[:-4], '*.pdb'))
            if sub_pdbs:
                logging.info("Picked up %d neglected PDBs: %s",
                             len(sub_pdbs), ' '.join(sub_pdbs))
                pdbs.extend(sub_pdbs)
            continue
        best_tmscore = -1
        best_pdb = None
        for rec in SeqIO.parse(sub_pdb_fname, 'fasta'):
            # Extract TM-score
            for token in rec.description.split():
                if token.startswith('TM-score'):
                    try:
                        this_tmscore = float(token.split('=', 1)[1])
                        if this_tmscore > best_tmscore:
                            best_tmscore = this_tmscore
                            best_pdb = rec.id
                    except:
                        logging.warn("PDB seq parsing issue: %s",
                                     rec.description)
                    finally:
                        break
        if best_pdb is None:
            logging.warn("Empty PDB alignment: " + sub_pdb_fname)
        else:
            logging.info("Best PDB of %s: %s", sub_pdb_fname, best_pdb)
            pdbs.append(best_pdb)

    records = tmalign.align_structs(pdbs,
                                    [seed for seed in map(str, sub_pdb_seqs)
                                    if isfile(seed) and not is_empty(seed)])
    SeqIO.write(records, task.target, 'fasta')
    # if not records:
    #     logging.info("Created empty PDB alignment %s", task.target)


def align_profiles(task, use_pdb=None):
    """Align several FASTA files with MAFFT. Clustal output.

    Cleans: [depends].cons.seq, [target].families.fa, [target].families.seq,
        [target].seq
    """
    seeds, singles = [], []
    # PDB alignment -- include as a seed profile if requested
    subalignments, pdb_seed = task.depends[:-1], str(task.depends[-1])
    if use_pdb and not is_empty(pdb_seed):
        seeds.append(pdb_seed)
    else:
        logging.info("Empty PDB alignment: %s", pdb_seed)

    # Get subfamily and subgroup consensus sequences/profiles
    for subaln in subalignments:
        aln = AlignIO.read(str(subaln), 'clustal')
        # with open(task.target, 'w+') as outfile:
        with open(ext(subaln, 'cons.seq'), 'w+') as outfile:
            outfile.write(">%s consensus\n" % basename(noext(subaln)))
            cons_seq = consensus.consensus(aln, trim_ends=False)
            if isdir(noext(subaln)):
                # Group profiles: include the subfamily consenses, too
                outfile.write(cons_seq + "\n")
                for record in aln:
                    outfile.write(">%s\n" % record.id)
                    outfile.write("%s\n" % record.seq)
            else:
                # Ungapped family consensus sequences
                outfile.write(cons_seq.replace('-', '') + "\n")

    # Merge the sequences and profiles
    for subconsseq in ext(subalignments, 'cons.seq'):
        if isdir(subconsseq[:-9]):
            # Group
            seeds.append(subconsseq)
        else:
            singles.append(subconsseq)
    # First, align/merge the single family consensus sequences
    famfa = ext(task.target, 'families.fa')
    allseq = ext(task.target, 'seq')
    assert singles or seeds, \
            'No .fasta files found to build %s' % task.target
    if singles:
        sh("cat %s > %s" % (' '.join(singles), famfa))
    if seeds:
        # Align the families with the groups
        sh("mafft --quiet --amino --globalgenafpair --maxiterate 1000 %s %s > %s"
           % (' '.join(['--seed '+s for s in seeds]), famfa, allseq))
        # XXX fast version
        # sh("mafft --quiet --amino --auto %s %s > %s"
        #    % (' '.join(['--seed '+s for s in seeds]), famfa, allseq))
    else:
        # No group profiles -- just align the families
        sh("mafft --quiet --amino --globalgenafpair --maxiterate 1000 %s > %s"
                % (famfa, allseq))
    # Convert FASTA to "pressed" (single-row) Clustal
    records = [rec for rec in SeqIO.parse(allseq, 'fasta')
            # Drop PDB-derived sequences
            # if ':' not in rec.id
            if 'TMalign' not in rec.description and
               'TM-score' not in rec.description
            ]
    records = list(alnutils.remove_empty_cols(records))
    if seeds:
        # MAFFT prefixes seed alignments with '_seed_' -- get rid of that
        for rec in records:
            if rec.id.startswith('_seed_'):
                rec.id = rec.id[6:]
    try:
        max_id_len = max(len(r.id) for r in records)
    except ValueError:
        # Common effup
        raise ValueError("Profile alignment failed for %s.\nInputs: %s"
                         % (task.target, ' '.join(map(str, task.depends))))

    with open(task.target, 'w+') as outfile:
        outfile.write('CLUSTAL X (-like) multiple sequence alignment\n\n')
        outfile.writelines(
                ['%s %s\n' % (rec.id.ljust(max_id_len), rec.seq)
                 for rec in records])


def aln2hmm(task):
    """Convert a Clustal alignment to an HMM profile.

    Cleans: .stk
    """
    stk = ext(task.depends[0], 'stk')
    SeqIO.convert(str(task.depends[0]), 'clustal', stk, 'stockholm')
    sh('hmmbuild %s %s' % (task.target, stk))


def cat_sub_consenses(task):
    """Concatenate the subfamily consensus sequences."""
    with open(task.target, 'w+') as outfile:
        for subaln in ext(task.depends, 'aln'):
            aln = AlignIO.read(str(subaln), 'clustal')
            outfile.write(">%s consensus\n" % noext(subaln))
            outfile.write(consensus.consensus(aln, trim_ends=False) + "\n")
            # Group profiles: include the subfamily consenses, too
            if isdir(noext(subaln)):
                with open(ext(subaln, 'fasta')) as subfam_file:
                    outfile.write(subfam_file.read())


def all_hmm(task):
    """Concatenate all HMM profiles into a database & press."""
    base = noext(task.depends[-1])
    # TODO - filter out *_all.hmm from `find` hits
    sh("cat %s.hmm `find %s/ -name '*.hmm' | grep -v '_all.hmm'` > %s"
       % (base, base, task.target))
    sh("hmmpress -f %s" % task.target)


# MAPGAPS actions

def mg_aln2cma(task, level=None):
    """Convert an alignment profile to CMA (or .tpl).

    Depends: .aln
    Cleans: .cons.cma, .cons_iron.cma
    """
    base = noext(task.target)
    name = basename(base)
    # Add consensus back to the subfamily-consensus seq set (.aln)
    # to produce a CMA (.cons.cma)
    aln = AlignIO.read(str(task.depends[0]), 'clustal')
    cons_rec = SeqRecord(Seq(consensus.consensus(aln, trim_ends=False)),
                         id=name, description=name + ' consensus')
    aln._records.insert(0, cons_rec)
    # Tidy up the CMA
    cmaln = biocma.ChainMultiAlignment(aln, level=level)
    biocma.write([cmaln], task.target, do_iron=True)

    # -------------
    ### OR (HMMer only):
    ### hmmemit consensus & reuse the .stk (done for .hmm) directly
    ###   see hmmalign --mapali option to include original .stk
    ###   .stk is no longer for 'clean'; original must be retained
    # stk = ext(task.target, 'stk')
    # cons_fa = ext(task.target, 'cons.fa')
    # sh('hmmalign --amino %s %s > %s'
    #         % (task.depends[1], cons_fa, stk))
    # SeqIO.convert(stk, 'stockholm', cons_fa, 'fasta')
    # sh("press < %s > %s" % (cons_fa, cons_seq))
    # sh("fa2cma %s > %s" % (cons_seq, base + '.cons.cma'))
    # -------------


def mg_cat_cma(task):
    """Concatenate subfamily MAPGAPS profiles.

    Depends: .cma of each subfamily, .tpl of this & sub-groups.
    """
    assert task.depends, 'No CMA files were given'
    sh('cat %s > %s' % (' '.join(map(str, task.depends)), task.target))


def mg_run_map(task):
    """Build/compile the complete set of MAPGAPS profiles."""
    sh("run_map %s" % noext(task.target))


# Tree

def treeify_subdirs(topdir):
    """Build a Newick string from the directory tree structure.

    Internal nodes = names of non-empty dirs
    External nodes = names of .fasta files
    """
    # Full paths of directories under topdir
    subdirs = filter(isdir, 
            [join(topdir, sd) for sd in sorted(os.listdir(topdir))])
    # Do internal nodes first, then the tips
    # Internal nodes: subtree Newick strings
    subtree_strs = [treeify_subdirs(sd) for sd in sorted(subdirs)]
    # Tips: string names (basename minus trailing .fasta)
    tip_names = [basename(fafname)[:-len('.fasta')]
            for fafname in glob(join(topdir, '*.fasta'))
            if noext(fafname) not in subdirs]
    tip_names.sort(key=str)
    return '(%s)%s' % (
            ','.join(subtree_strs + tip_names),
            basename(topdir.rstrip('/')))



# === Scan ============================================
# ENH - use the built HMM/MAPGAP & Newick tree to scan&classify seqs

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

    # Write a table of just the best hit for each sequence
    if args.table:
        write_table(scanhits)
    else:
        write_summary(scanhits)
    if args.seqout:
        write_seqout(scanhits, args.target, args.seqout)
    if args.seqsets:
        write_seqsets(scanhits, args.target)


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



# === Add ============================================
# Leave-one-out validation and refinement of each profile

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



# === Refine ============================================
# Leave-one-out validation and refinement of each profile

def cmd_refine(args):
    assert isdir(args.basedir)

    # Defaults: locate args.unclassified or set it
    if args.unclassified:
        unclass_fa = args.unclassified
    else:
        topname = basename(abspath(args.basedir).rstrip('/'))
        unclass_fa = join(args.basedir, topname + '-Unique.fasta')
    if not isfile(unclass_fa):
        raise ValueError("Can't find the unclassified sequence set "
                            "(%s does not exist)" % unclass_fa)

    all_fa = glob(join(args.basedir, '*.fasta'))
    all_fa.remove(unclass_fa)
    if args.exclude:
        if args.exclude in all_fa:
            all_fa.remove(args.exclude)
        elif args.exclude + '.fasta' in all_fa:
            all_fa.remove(args.exclude + '.fasta')
        else:
            logging.warn("Didn't find excluded profile %s in %s",
                         args.exclude, args.basedir)
    # Prerequisites: HMM profiles for each sequence set must be built
    for fa in all_fa:
        if not isfile(ext(fa, 'hmm')):
            raise ValueError("Unbuilt HMM for sequence set: %s" % fa)

    for this_fa in all_fa:
        this_hmm = ext(this_fa, 'hmm')
        logging.info("Refining %s...", this_fa)

        # Scan the seed sequences in `this_fa` w/ `this_hmm`
        # -> {seq_id: top_dom_eval}
        seq_evals = get_hmm_evalues(this_hmm, this_fa)

        # Starting with the lowest-scoring seed, do leave-one-out.
        # Leave at least 2 sequences in the profile, otherwise MAFFT will fail
        # (and a set of only 2 seqs should be merged into Unclassified anyway)
        for seq_id, evalue in sorted(seq_evals.items()[:-2],
                                    key=lambda kv: kv[1]):
            # Copy all seqs in this_fa, except seq_id, to a temp file
            with tempfile.NamedTemporaryFile(suffix='.fa') as tempseqfile:
                for rec in SeqIO.parse(this_fa, 'fasta'):
                    if rec.id == seq_id:
                        # Put this seq in a separate file
                        SeqIO.write(rec, tempseqfile.name + '.seed', 'fasta')
                    else:
                        SeqIO.write(rec, tempseqfile, 'fasta')
                tempseqfile.flush()
                # Align the "lean" profile w/ MAFFT (as usual)
                sh("mafft --quiet --amino --reorder --maxiterate 1000 "
                   "--genafpair --ep 0.123 %s > %s.seq"
                   % (tempseqfile.name, tempseqfile.name))
                # Build the "lean" HMM w/ hmmbuild
                SeqIO.convert(tempseqfile.name + '.seq', 'fasta',
                              tempseqfile.name + '.stk', 'stockholm')
                sh('hmmbuild %s %s > /dev/null'
                   % (tempseqfile.name + '.hmm', tempseqfile.name + '.stk'))
                os.remove(tempseqfile.name + '.seq')
                os.remove(tempseqfile.name + '.stk')

                # Scan the dropped seed seq w/ the "lean" HMM
                seed_eval = get_hmm_evalues(tempseqfile.name + '.hmm',
                                            tempseqfile.name + '.seed'
                                           )[seq_id]
                # Scan the unclassified seqs w/ the "lean" HMM
                # Keep the top evalue only; that's our threshold
                db_evals = get_hmm_evalues(tempseqfile.name + '.hmm',
                                           unclass_fa)
                try:
                    out_eval = max(db_evals.values())
                except ValueError:
                    # No HMM hits!
                    logging.warn("No unclassified hits for %s!", this_fa)
                    break

                if out_eval <= seed_eval:
                    if args.dry_run:
                        logging.info("* drop %s (evalue %s vs. %s)",
                                     seq_id, seed_eval, out_eval)
                    else:
                        # TODO:
                        # - append the dropped seq to unclass_fa
                        # - replace this_fa w/ the "lean" profile
                        # - Logging: print the name of the dropped seq
                        logging.info(
                            "Dropped %s from %s into %s (evalue %s vs. %s)",
                            seq_id, this_fa, unclass_fa, seed_eval, out_eval)
                else:
                    # Since the seeds are sorted by evalue, if this sequence is a
                    # keeper, then so are the rest. On to the next profile!
                    if args.dry_run:
                        logging.info("* keep %s (evalue %s vs. %s)",
                                     seq_id, seed_eval, out_eval)
                    break


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


# === Update FASTA ============================================
# Overwrite the original FASTA files with ungapped .aln sequences

def cmd_update_fasta(args):
    if isdir(args.basedir):
        for task in find_update_tasks(args.basedir):
            task.build()
    elif isfile(args.basedir) and args.basedir.endswith('.aln'):
        task = one_update_task(args.basedir)
        if task is not None:
            task.build()
    else:
        raise RuntimeError("Target must be a directory or .aln file: %s"
                         % args.basedir)


def find_update_tasks(topdir):
    """Locate and update all .fasta sequences."""
    # if <base>.aln and <base>.fasta exist (same base)
    # and .aln is newer than .fasta (handled by Task)
    # [ENH: and same # seqs?]
    # then update_fasta()
    for dirpath, _subdirs, fnames in os.walk(topdir):
        todo_bases = [join(dirpath, fn[:-4])
                      for fn in fnames
                      if fn.endswith('.aln') and fn[:-4] + '.fasta' in fnames]
        for base in todo_bases:
            if not is_same_aln_fasta(base + '.aln', base + '.fasta'):
                yield Task(base + '.fasta', 
                    action=update_fasta,
                    depends=base + '.aln')


def one_update_task(aln_fname):
    fasta_fname = aln_fname[:-4] + '.fasta'
    if not isfile(fasta_fname):
        raise RuntimeError("No corresponding .fasta file for %s: need %s"
                           % (aln_fname, fasta_fname))
    if is_same_aln_fasta(aln_fname, fasta_fname):
        logging.info("%s already matches %s." % (fasta_fname, aln_fname))
    else:
        return Task(fasta_fname, 
                    action=update_fasta,
                    depends=aln_fname)


def update_fasta(task):
    """Overwrite the .fasta with ungapped unique .aln sequences."""
    # ENH: warn about:
    # - Sequences containing 'X' (log the number of X's)
    # - Short sequences (e.g. if seq_length < .5 * alignment_width)

    def ungap_and_unique(records):
        """Remove gap characters, drop duplicate SeqRecords, and sort by length.

        Equivalent to:
        seqioconvert clustal fasta < t.aln | purge.py -c | press.py -u > t.fasta
        """
        seen = set()
        counter = 0
        for rec in sorted(records,
                key=lambda r: len(r.seq) - r.seq.count('-') - r.seq.count('.'),
                reverse=True):
            # Remove gaps ('-')
            thisseq = str(rec.seq).replace('-', '').replace('.', '').upper()
            # Drop if this sequence is completely contained by another
            for longseq in seen:
                if thisseq in longseq:
                    counter += 1
                    break
            else:
                seen.add(thisseq)
                # Apply the ungapping from above
                rec.seq._data = thisseq
                yield rec

        # I/O
        msg = "Updated %s" % task.target
        if counter:
            msg += (" (purged %d sequence%s)"
                    % (counter, '' if counter == 1 else 's'))
        logging.info(msg)

    records = ungap_and_unique(SeqIO.parse(str(task.depends[0]), 'clustal'))
    SeqIO.write(records, task.target, 'fasta')


def is_same_aln_fasta(aln_fn, fasta_fn):
    """Compare ungapped seqs in aln and fasta; are they all the same?"""
    def rec2tuple(rec, do_ungap):
        if do_ungap:
            strseq = str(rec.seq).replace('-', '').replace('.', '').upper()
        else:
            strseq = str(rec.seq).upper()
        return (rec.id, strseq)

    aln_hash = set([rec2tuple(rec, True)
                    for rec in SeqIO.parse(aln_fn, 'clustal')])
    fa_hash = set([rec2tuple(rec, False)
                   for rec in SeqIO.parse(fasta_fn, 'fasta')])
    return aln_hash == fa_hash


# === Main ============================================

if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(
            description=__doc__,
            epilog="Contact Eric <etal@uga.edu> for help.")
    # Global options
    AP.add_argument('--quiet',
            action='store_true',
            help="Don't print status messages, only warnings and errors.")
    AP_subparsers = AP.add_subparsers(
            help='sub-commands (use with -h for more info)')
    # Subcommand: build
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
    # Subcommand: update-fasta
    P_update = AP_subparsers.add_parser('update-fasta', 
            help="""Replace original FASTA sequence sets with the ungapped
            sequences from the corresponding alignment (.aln) files, sorted by
            decreasing length.""")
    P_update.add_argument('basedir',
            help='Top of the directory tree where the sequences sets are.')
    P_update.set_defaults(func=cmd_update_fasta)

    args = AP.parse_args()
    # Handle global options here
    if not args.quiet:
        logging.basicConfig(level=logging.INFO,
                format="@%(lineno)s | %(message)s")
    # Pass the rest along to the sub-command implementation
    args.func(args)

