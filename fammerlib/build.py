"Build alignments and profiles in a directory tree of sequence sets."

import itertools
import logging
import os
import subprocess
from glob import glob
from os.path import basename, isdir, isfile, join

from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biocma import biocma
from biofrills import consensus, alnutils

from ._share import write_fasta
from .tasks import Task, ext, noext, sh, which, is_empty
from . import tmalign


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
                action=align_fasta_mafft,
                depends=subfa,
                cleans=ext(subfa, 'seq')))  # mafft
                # cleans=[subfa + '.1.fas', subfa + '.2.fas']))  # prank
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

def align_fasta_mafft(task):
    """Align a FASTA file with MAFFT. Clustal output."""
    seq = ext(task.depends[0], 'seq')
    sh("mafft --quiet --amino --reorder --maxiterate 1000 "
       "--genafpair --ep 0.123 %s > %s"
       % (task.depends[0], seq))
    # Convert FASTA to "pressed" (single-row) Clustal
    records = list(SeqIO.parse(seq, 'fasta'))
    # Check for 'X' characters in the sequences -- these cause problems
    for rec in records:
        if 'X' in str(rec.seq):
            logging.warn('Sequence %r contains unknown residue X', rec.id)
    max_id_len = max(len(r.id) for r in records)
    with open(task.target, 'w+') as outfile:
        outfile.write('CLUSTAL X (-like) multiple sequence alignment\n\n')
        outfile.writelines(
                ['%s %s\n' % (rec.id.ljust(max_id_len), rec.seq)
                 for rec in records])


def align_fasta_prank(task):
    """Align a FASTA file with PRANK. Clustal output.

    Cleans: [input].fasta.{1,2}.fas
    """
    seq = task.depends[0] + '.2.fas'
    sh("prank -d=%s -o=%s -twice -quiet" % (task.depends[0], task.depends[0]))
    # Convert FASTA to "pressed" (single-row) Clustal
    records = list(SeqIO.parse(seq, 'fasta'))
    # Check for 'X' characters in the sequences -- these cause problems
    for rec in records:
        if 'X' in str(rec.seq):
            logging.warn('Sequence %r contains unknown residue X', rec.id)
    max_id_len = max(len(r.id) for r in records)
    with open(task.target, 'w+') as outfile:
        outfile.write('CLUSTAL X (-like) multiple sequence alignment\n\n')
        outfile.writelines(
                ['%s %s\n' % (rec.id.ljust(max_id_len), rec.seq)
                 for rec in records])


def align_pdbs(task, sub_pdb_seqs=(), use_pdb=None):
    """Create a structure-based sequence alignment from PDB files.

    Inputs are PDB files and FASTA alignments (of previously aligned PDBs).
    """
    if not use_pdb or not task.depends:
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
                    except Exception:
                        logging.warn("PDB seq parsing issue: %s",
                                     rec.description)
                    finally:
                        break
        if best_pdb is None:
            logging.warn("Empty PDB alignment: " + sub_pdb_fname)
        else:
            logging.info("Best PDB of %s: %s", sub_pdb_fname, best_pdb)
            pdbs.append(best_pdb)

    pdbseedfnames = [seed for seed in map(str, sub_pdb_seqs)
                     if isfile(seed) and not is_empty(seed)]

    try:
        mustang_tmpfname = '_tmp_mustang.afasta'
        if len(pdbs) > 1 and which(['mustang']):
            # Align PDBs with MUSTANG.
            subprocess.check_call(['mustang',
                                    '-o', '_tmp_mustang',
                                    '-F', 'fasta',
                                    '-s', 'OFF',
                                    '-i'] + pdbs)
            pdbseedfnames.append(mustang_tmpfname)

        # This is where the magic happens.
        records = tmalign.align_structs(pdbs, pdbseedfnames)

    finally:
        if isfile(mustang_tmpfname):
            os.remove(mustang_tmpfname)

    # SeqIO.write(records, task.target, 'fasta')
    write_fasta(records, task.target)
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
            cons_seq = consensus.consensus(aln, trim_ends=False,
                                           gap_threshold=0.6)
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
               'TM-score' not in rec.description and
               not rec.id.endswith('.pdb')
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
            outfile.write(consensus.consensus(aln, trim_ends=False,
                                              gap_threshold=0.6) + "\n")
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
    cons_rec = SeqRecord(Seq(consensus.consensus(aln, trim_ends=False,
                                                 gap_threshold=0.6)),
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


