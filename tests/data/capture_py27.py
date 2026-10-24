#!/usr/bin/env python2.7
"""Capture reference outputs of the Python 2 fammer code base.

This script must run under Python 2.7 with Biopython 1.76, against the
unmodified fammer sources of FAMMER_SOURCES_COMMIT below, with MAFFT, HMMER,
FastTree and TM-align on PATH; see README.md beside it for the environment.
biofrills and biocma are taken from their sibling Git checkouts at the pinned
Python 2 commits below, not from their working trees.  The capture puts
tool-shims/ first on PATH; its mafft wrapper replaces --globalgenafpair,
which current MAFFT rejects, with --globalpair.

It rewrites two directories beside this file:

* py27/ -- what the Python 2 code computed:
  - library/*.json: pure-Python functions applied to synthetic and real
    inputs, and to the tool outputs below;
  - library/mapgaps/: the .cma and .tpl files written by build.mg_aln2cma;
  - e2e/*.json: command-line runs of fammer.py and tmalign.py on a scratch
    copy of profiles/SomeKinases, with stdout, stderr, exit status and the
    files each run wrote or changed;
  - e2e/build/: the alignment-level files written by `fammer build`.
* tool-outputs/ -- raw output of the external programs, for testing the
  parsers without the programs installed.

Every function call records either {"value": ...} or {"error": {"type": ...,
"message": ...}}, so failures of the original code are part of the reference.
Floats are written with Python 2.7's repr, the shortest decimal string that
reads back to the identical double.

Usage: see "Regenerating" in README.md.  The capture does not run from the
repository root; the comment on FAMMER_SOURCES_COMMIT below says why.
"""
from __future__ import print_function

import argparse
import hashlib
import io
import json
import logging
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(HERE, 'py27')
TOOL_OUT = os.path.join(HERE, 'tool-outputs')
SHIMS = os.path.join(HERE, 'tool-shims')

# The last Python 2 commits of the two upstream libraries.  biofrills
# 36684bb is its last commit before the Python 3 port; biocma eac0c57 is the
# commit biofrills' own reference capture uses.
BIOFRILLS_COMMIT = '36684bb6c7632f96215e8b2b4ebc86640f331bcd'
BIOCMA_COMMIT = 'eac0c57eb83a9498e53ccdeb9cbc3fe21a5826a7'

# fammer's code under test: the commit of SOURCES whose outputs the references
# record.  Any later change to these sources, even one that preserves behavior,
# moves the source line numbers that the recorded log lines carry, so the
# capture refuses any other sources; README.md describes running it from a
# worktree.
FAMMER_SOURCES_COMMIT = '5bcf8ebcfcb4242d81df896119e97741a1bb9f08'
SOURCES = ['fammerlib', 'fammer.py', 'tmalign.py', 'setup.py']

# Input data copied from the repository; checked and recorded like SOURCES.
INPUTS = ['profiles', 'test']

TOOLS = ['mafft', 'hmmbuild', 'hmmsearch', 'hmmpress', 'fasttree', 'TMalign']

# TM-align pairs whose raw output is kept, as paths relative to the repository.
TMALIGN_PAIRS = [
    ('profiles/ROPK/1ATP.pdb', 'profiles/ROPK/1CSN.pdb'),
    ('profiles/ROPK/2W1Z.pdb', 'profiles/ROPK/3Q60.pdb'),
    ('profiles/SomeKinases/CMGC/1JST.pdb', 'profiles/SomeKinases/CMGC/2ERK.pdb'),
]


# ---------------------------------------------------------------------------
# Helpers

def outcome(func, *args, **kwargs):
    """Call func; return {"value": result} or {"error": {...}}."""
    try:
        return OrderedDict([('value', func(*args, **kwargs))])
    except Exception as exc:
        return OrderedDict([('error', OrderedDict([
            ('type', type(exc).__name__),
            ('message', str(exc)),
        ]))])


def records_table(records):
    """SeqRecords as [id, description, sequence] rows."""
    return [[rec.id, rec.description, str(rec.seq)] for rec in records]


def make_records(rows):
    """[id, sequence] or [id, description, sequence] rows as SeqRecords."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    records = []
    for row in rows:
        if len(row) == 2:
            records.append(SeqRecord(Seq(row[1]), id=row[0]))
        else:
            records.append(SeqRecord(Seq(row[2]), id=row[0],
                                     description=row[1]))
    return records


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def read_bytes(path):
    with open(path, 'rb') as handle:
        return handle.read()


def write_bytes(path, data):
    if not os.path.isdir(os.path.dirname(path)):
        os.makedirs(os.path.dirname(path))
    with open(path, 'wb') as handle:
        handle.write(data)


def normalize_strings(data):
    """data with NORMALIZE applied to every string in it."""
    if isinstance(data, basestring):  # noqa: F821 (Python 2 builtin)
        return NORMALIZE(data)
    if isinstance(data, dict):
        return OrderedDict((key, normalize_strings(value))
                           for key, value in data.items())
    if isinstance(data, (list, tuple)):
        return [normalize_strings(value) for value in data]
    return data


def write_json(relpath, data, root=OUT):
    """Write data under root as indented ASCII JSON, with machine-specific
    paths replaced by labels."""
    text = json.dumps(normalize_strings(data), indent=1,
                      separators=(',', ': '), ensure_ascii=True)
    write_bytes(os.path.join(root, relpath), text.encode('ascii') + b'\n')


class LogRecorder(logging.Handler):
    """Collect log records as [level, module, message] rows."""

    def __init__(self):
        logging.Handler.__init__(self)
        self.rows = []

    def emit(self, record):
        self.rows.append([record.levelname, record.module,
                          record.getMessage()])


def with_log(func, *args, **kwargs):
    """outcome() of a call, plus the INFO-and-above messages it logged."""
    recorder = LogRecorder()
    root = logging.getLogger()
    saved_level = root.level
    root.setLevel(logging.INFO)
    root.addHandler(recorder)
    try:
        result = outcome(func, *args, **kwargs)
    finally:
        root.removeHandler(recorder)
        root.setLevel(saved_level)
    result['log'] = recorder.rows
    return result


class Normalizer(object):
    """Replace machine-specific paths in recorded text."""

    def __init__(self):
        self.pairs = []

    def add(self, path, label):
        for variant in (os.path.realpath(path), path):
            if variant not in [p for p, _ in self.pairs]:
                self.pairs.append((variant, label))
        # Longest first, so a prefix does not shadow a longer path.
        self.pairs.sort(key=lambda pair: -len(pair[0]))

    def __call__(self, text):
        for path, label in self.pairs:
            text = text.replace(path, label)
        # Random names of temporary files and directories.
        return re.sub(r'<tmpdir>/tmp[\w.-]+', '<tmpfile>', text)


NORMALIZE = Normalizer()


def git_archive(repo, commit, path, destdir):
    """Extract path at commit from the Git checkout repo into destdir."""
    archive = subprocess.Popen(['git', '-C', repo, 'archive', commit, path],
                               stdout=subprocess.PIPE)
    subprocess.check_call(['tar', '-x', '-C', destdir], stdin=archive.stdout)
    archive.stdout.close()
    assert archive.wait() == 0, "git archive of %s failed" % repo


def copy_tracked(relpath, destdir):
    """Copy the Git-tracked files under relpath (repository-relative) into
    destdir, keeping their paths relative to relpath's parent directory."""
    names = subprocess.check_output(['git', 'ls-files', '--', relpath],
                                    cwd=REPO).decode('utf-8').splitlines()
    parent = os.path.dirname(relpath)
    for name in names:
        dest = os.path.join(destdir, os.path.relpath(name, parent))
        if not os.path.isdir(os.path.dirname(dest)):
            os.makedirs(os.path.dirname(dest))
        shutil.copy(os.path.join(REPO, name), dest)


def snapshot(topdir):
    """{relative path: sha256} of every file under topdir."""
    result = {}
    for dirpath, _dirs, fnames in os.walk(topdir):
        for fname in fnames:
            path = os.path.join(dirpath, fname)
            result[os.path.relpath(path, topdir)] = sha256(read_bytes(path))
    return result


def tool_version(argv, pattern):
    """First line of a program's output matching pattern."""
    proc = subprocess.Popen(argv, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    output = proc.communicate()[0].decode('utf-8', 'replace')
    for line in output.splitlines():
        if re.search(pattern, line):
            return line.strip()
    return None


# ---------------------------------------------------------------------------
# Library: tasks

def capture_tasks():
    from fammerlib import tasks
    result = OrderedDict()

    paths = ['fam.fasta', 'dir/fam.fasta', 'dir.v2/fam', 'noext', 'a.b.c',
             'tree/.hidden', ['x.fasta', 'y/z.aln']]
    result['ext'] = [[p, outcome(tasks.ext, p, 'aln')] for p in paths]
    result['noext'] = [[p, outcome(tasks.noext, p)] for p in paths]
    task = tasks.Task('group/fam.aln')
    result['ext_of_task'] = outcome(tasks.ext, task, 'hmm')
    result['noext_of_task'] = outcome(tasks.noext, task)

    workdir = tempfile.mkdtemp(prefix='fammer-tasks-')
    NORMALIZE.add(workdir, '<taskdir>')
    try:
        sizes = OrderedDict()
        for size in (0, 1, 2, 10):
            path = os.path.join(workdir, 'size%d' % size)
            write_bytes(path, b'x' * size)
            sizes[str(size)] = outcome(tasks.is_empty, path)
        result['is_empty_by_size'] = sizes
        result['build_scenarios'] = capture_task_scenarios(workdir)
    finally:
        shutil.rmtree(workdir)

    # sh(): success, failure without output, failure with output.
    result['sh'] = [
        [cmd, outcome(tasks.sh, cmd)]
        for cmd in ('true', 'exit 3', 'echo to-stdout; echo to-stderr >&2; exit 2')]
    result['which_missing'] = outcome(tasks.which,
                                      ['fammer-no-such-program'])
    return result


def capture_task_scenarios(workdir):
    """Task.build/is_up_to_date/clean over files with controlled mtimes.

    Every initial file gets an mtime in 2001; a rebuilt target gets the real
    current time, so it is newer than any initial file.
    """
    from fammerlib.tasks import Task
    t_old, t_mid, t_new = 1000000000, 1000000100, 1000000200

    def path(name):
        return os.path.join(workdir, name)

    def make(name, mtime):
        write_bytes(path(name), b'content\n')
        os.utime(path(name), (mtime, mtime))

    calls = []

    def action(task):
        calls.append(os.path.basename(task.target))
        write_bytes(task.target, b'built\n')

    def failing_action(task):
        calls.append(os.path.basename(task.target))
        raise ValueError('action failed')

    def reset():
        for fname in os.listdir(workdir):
            os.remove(path(fname))
        del calls[:]

    scenarios = OrderedDict()

    def run(name, setup, top):
        reset()
        setup()
        up_to_date = outcome(top.is_up_to_date)
        built = with_log(top.build)
        for row in built['log']:
            row[2] = row[2].replace(workdir + os.sep, '')
        scenarios[name] = OrderedDict([
            ('is_up_to_date_before', up_to_date),
            ('build', built),
            ('actions_called', list(calls)),
            ('is_up_to_date_after', outcome(top.is_up_to_date)),
        ])

    run('target_missing',
        lambda: make('dep.in', t_old),
        Task(path('out'), action, depends=path('dep.in')))
    run('target_newer_than_dep',
        lambda: (make('dep.in', t_old), make('out', t_new)),
        Task(path('out'), action, depends=path('dep.in')))
    run('target_same_mtime_as_dep',
        lambda: (make('dep.in', t_mid), make('out', t_mid)),
        Task(path('out'), action, depends=path('dep.in')))
    run('target_older_than_dep',
        lambda: (make('dep.in', t_new), make('out', t_old)),
        Task(path('out'), action, depends=path('dep.in')))
    run('stale_middle_task',
        lambda: (make('dep.in', t_mid), make('mid', t_old),
                 make('top', t_new)),
        Task(path('top'), action,
             depends=Task(path('mid'), action, depends=path('dep.in'))))
    run('missing_file_dependency_with_target',
        lambda: make('out', t_new),
        Task(path('out'), action, depends=path('absent.in')))
    run('missing_file_dependency_without_target',
        lambda: None,
        Task(path('out'), action, depends=path('absent.in')))
    run('failing_action',
        lambda: make('dep.in', t_old),
        Task(path('out'), failing_action, depends=path('dep.in')))

    # clean(): removes listed files that exist, recursively through tasks.
    reset()
    for name in ('a.tmp', 'b.tmp', 'c.tmp', 'out', 'mid'):
        make(name, t_old)
    top = Task(path('out'), cleans=[path('a.tmp'), path('absent.tmp')],
               depends=[Task(path('mid'), cleans=path('b.tmp')),
                        path('c.tmp')])
    cleaned = with_log(top.clean)
    for row in cleaned['log']:
        row[2] = row[2].replace(workdir + os.sep, '')
    scenarios['clean'] = OrderedDict([
        ('clean', cleaned),
        ('remaining_files', sorted(os.listdir(workdir))),
    ])
    reset()
    return scenarios


# ---------------------------------------------------------------------------
# Library: build (task graph and Newick tree, no external programs)

def describe_dep(dep):
    from fammerlib.tasks import Task
    if isinstance(dep, Task):
        return OrderedDict([('task', dep.target)])
    return OrderedDict([('file', dep)])


def task_graph(result):
    """Every task reachable from a build.Result, keyed by target."""
    from fammerlib.tasks import Task
    tasks = {}
    conflicts = []

    def visit(task):
        if not isinstance(task, Task):
            return
        entry = OrderedDict([
            ('action', task.action.__name__),
            ('kwargs', OrderedDict(sorted(task.kwargs.items()))),
            ('depends', [describe_dep(dep) for dep in task.depends]),
            ('cleans', list(task.cleans)),
        ])
        if task.target in tasks:
            if tasks[task.target] != entry:
                conflicts.append(task.target)
            return
        tasks[task.target] = entry
        for dep in task.depends:
            visit(dep)

    roots = OrderedDict()
    for attr in ('aln', 'pdbseq', 'hmm', 'cma', 'tpl'):
        task = getattr(result, attr)
        roots[attr] = task.target if isinstance(task, Task) else task
        visit(task)
    return OrderedDict([
        ('result', roots),
        ('tasks', OrderedDict(sorted(tasks.items()))),
        ('conflicting_targets', sorted(set(conflicts))),
    ])


# Group FASTA files placed beside the group directories of SomeKinases, a
# layout the bundled trees lack.  taskify_subdirs and treeify_subdirs skip a
# FASTA file whose name matches a sibling directory; this is the code path
# that a Python 3 port breaks if `subdirs` becomes a one-shot iterator.
GROUP_FASTAS = ['CAMK.fasta', 'CMGC.fasta', 'STE.fasta']


def capture_build_graph(tree, group_fastas=False):
    from fammerlib import build
    workdir = tempfile.mkdtemp(prefix='fammer-graph-')
    cwd = os.getcwd()
    try:
        copy_tracked(os.path.join('profiles', tree), workdir)
        if group_fastas:
            for fname in GROUP_FASTAS:
                write_bytes(os.path.join(workdir, tree, fname),
                            b'>placeholder\nACDEFGHIK\n')
        os.chdir(workdir)
        result = OrderedDict()
        result['tree'] = 'profiles/' + tree
        if group_fastas:
            result['added_files'] = [tree + '/' + fname
                                     for fname in GROUP_FASTAS]
        result['treeify_subdirs'] = outcome(build.treeify_subdirs, tree)
        graphs = []
        for hmmer in (False, True):
            for mapgaps in (False, True):
                for use_pdb in (False, True):
                    graphs.append(OrderedDict([
                        ('kwargs', OrderedDict([('hmmer', hmmer),
                                                ('mapgaps', mapgaps),
                                                ('use_pdb', use_pdb),
                                                ('level', 0)])),
                        ('graph', outcome(lambda: task_graph(
                            build.taskify_subdirs(tree, hmmer, mapgaps,
                                                  use_pdb, 0)))),
                    ]))
        result['taskify_subdirs'] = graphs
        return result
    finally:
        os.chdir(cwd)
        shutil.rmtree(workdir)


# ---------------------------------------------------------------------------
# Library: build, MAPGAPS profiles (biocma; the MAPGAPS programs are not run)

# One case per level of the SomeKinases profile tree, with the target and
# level that taskify_subdirs assigns: a family profile, a group template and
# the top-level template.  The alignments are those of the end-to-end build.
MAPGAPS_CASES = [
    ('SomeKinases/CAMK/CDPK.aln', 'SomeKinases/CAMK/CDPK.cma', 2),
    ('SomeKinases/CAMK.aln', 'SomeKinases/CAMK.tpl', 1),
    ('SomeKinases.aln', 'SomeKinases.tpl', 0),
]


def capture_mapgaps(build_root):
    """build.mg_aln2cma on alignments that the end-to-end build wrote.

    capture_e2e keeps those alignments under build_root, so this runs after
    it.  Each written file is kept under library/mapgaps/ at its path in the
    build tree.  The consensus sequence that mg_aln2cma prepends is recorded
    too, recomputed here with the arguments of the consensus call inside
    mg_aln2cma (trim_ends=False, gap_threshold=0.6), so that biocma's writer
    can be checked apart from the biofrills consensus.
    """
    from Bio import AlignIO
    from biofrills import consensus
    from fammerlib import build
    from fammerlib.tasks import Task
    workdir = tempfile.mkdtemp(prefix='fammer-mapgaps-')
    NORMALIZE.add(workdir, '<mapgapsdir>')
    cases = []
    try:
        for aln_name, target_name, level in MAPGAPS_CASES:
            aln_path = os.path.join(build_root, aln_name)
            target = os.path.join(workdir, target_name)
            if not os.path.isdir(os.path.dirname(target)):
                os.makedirs(os.path.dirname(target))
            task = Task(target, action=build.mg_aln2cma, depends=aln_path,
                        kwargs={'level': level})
            before = snapshot(workdir)
            result = outcome(build.mg_aln2cma, task, **task.kwargs)
            created = sorted(set(snapshot(workdir)) - set(before))
            if os.path.isfile(target):
                write_bytes(os.path.join(OUT, 'library', 'mapgaps', target_name),
                            read_bytes(target))
            aln = AlignIO.read(aln_path, 'clustal')
            cases.append(OrderedDict([
                ('alignment', os.path.relpath(aln_path, OUT)),
                ('target', target_name),
                ('level', level),
                ('consensus', consensus.consensus(aln, trim_ends=False,
                                                  gap_threshold=0.6)),
                ('mg_aln2cma', result),
                ('created', created),
            ]))
    finally:
        shutil.rmtree(workdir)
    return OrderedDict([
        ('cases', cases),
        ('kept_files', 'library/mapgaps/: each target at its path in the '
                       'build tree'),
    ])


# ---------------------------------------------------------------------------
# Library: _share, scan, add

TBLOUT_DOCSTRING_EXAMPLE = """\
#                                                                  --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name           accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#   ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
CDPK-consensus          -          Testfam              -           1.9e-153  508.4  25.4  2.1e-153  508.2  17.6   1.0   1   0   0   1   1   1   1 -
gi|254575025|pdb|3HX4|A -          Testfam              -           1.4e-150  498.9   5.7  1.6e-150  498.7   4.0   1.0   1   0   0   1   1   1   1 {<Apicomplexa(e)>}Chain A, Crystal Structure Of Cdpk1 Of Toxoplasma Gondii, Tgme49_1014 Presence Of Calcium
TAIR|locus:2011201      -          Testfam              -           2.4e-137  455.2   0.0  2.7e-137  455.1   0.0   1.0   1   0   0   1   1   1   1 symbol:CDPK1 species:3702 "Arabidopsis thaliana"  AT1G18890
"""

# The same sequence hit by several profiles: lower, higher and equal domain
# E-values, in both orders.
TBLOUT_SYNTHETIC = """\
# synthetic
seqA - ProfOne - 1e-50 100.0 0.0 2e-50 99.0 0.0 1.0 1 0 0 1 1 1 1 -
seqA - ProfTwo - 1e-60 110.0 0.0 3e-60 109.0 0.0 1.0 1 0 0 1 1 1 1 -
seqA - ProfThree - 1e-55 105.0 0.0 5e-50 98.0 0.0 1.0 1 0 0 1 1 1 1 -
seqB - ProfTwo - 1e-10 30.0 0.0 1e-10 30.0 0.0 1.0 1 0 0 1 1 1 1 two words
seqB - ProfOne - 1e-10 30.0 0.0 1e-10 30.0 0.0 1.0 1 0 0 1 1 1 1 -
seqC - ProfOne - 0.5 1.0 0.0 0.5 1.0 0.0 1.0 1 0 0 1 1 1 1 -
"""


def hits_table(hits):
    """{seq: (profile, evalue)} as sorted [seq, profile, evalue] rows."""
    return [[seq, prof, evalue] for seq, (prof, evalue) in sorted(hits.items())]


def capture_stdout(func, *args, **kwargs):
    """outcome() of a call, with what it printed as 'stdout'."""
    saved = sys.stdout
    sys.stdout = io.BytesIO()
    try:
        result = outcome(func, *args, **kwargs)
        result['stdout'] = sys.stdout.getvalue().decode('utf-8')
    finally:
        sys.stdout = saved
    return result


def capture_share(tbl_path):
    from Bio import AlignIO, SeqIO
    from fammerlib import _share, add, scan
    result = OrderedDict()

    sources = OrderedDict([
        ('docstring_example', TBLOUT_DOCSTRING_EXAMPLE),
        ('synthetic', TBLOUT_SYNTHETIC),
        ('tool-outputs/' + os.path.relpath(tbl_path, TOOL_OUT),
         read_bytes(tbl_path).decode('utf-8')),
    ])
    parsed = OrderedDict()
    for name, text in sources.items():
        lines = text.encode('utf-8').splitlines(True)
        parsed[name] = outcome(lambda: hits_table(
            _share.parse_scanned_hits(lines)))
    result['inputs'] = OrderedDict([('docstring_example',
                                     TBLOUT_DOCSTRING_EXAMPLE),
                                    ('synthetic', TBLOUT_SYNTHETIC)])
    result['parse_scanned_hits'] = parsed

    real_hits = _share.parse_scanned_hits(read_bytes(tbl_path).splitlines(True))
    synthetic_hits = _share.parse_scanned_hits(
        TBLOUT_SYNTHETIC.splitlines(True))

    result['scan.write_summary'] = OrderedDict(
        (name, capture_stdout(scan.write_summary, hits))
        for name, hits in (('synthetic', synthetic_hits),
                           ('real', real_hits)))
    result['scan.write_table'] = OrderedDict(
        (name, capture_stdout(scan.write_table, hits))
        for name, hits in (('synthetic', synthetic_hits),
                           ('real', real_hits)))

    profiles = sorted(set(prof for prof, _e in real_hits.values()))
    result['add.filter_profile'] = OrderedDict(
        (prof, outcome(lambda prof=prof: sorted(
            add.filter_profile(real_hits, prof))))
        for prof in profiles + ['NoSuchProfile'])
    evalues = sorted(set(e for _p, e in real_hits.values()))
    thresholds = [0.0, evalues[len(evalues) // 2], evalues[-1], 1.0]
    result['add.filter_evalue'] = [
        [threshold, outcome(lambda threshold=threshold: list(
            add.filter_evalue(sorted(add.filter_profile(
                real_hits, profiles[0])), threshold)))]
        for threshold in thresholds]
    result['add.filter_evalue_profile'] = profiles[0]

    # write_fasta: a FASTA-parsed record (description includes the ID), a
    # Clustal-parsed record, and records built without a description.
    fasta_rec = next(SeqIO.parse(io.StringIO(
        u'>seq1 some description\nACDEFGHIKLMNPQRSTVWY\n'), 'fasta'))
    clustal_rec = AlignIO.read(
        os.path.join(REPO, 'profiles', 'ROPK', 'PKL.aln'), 'clustal')[0]
    built = make_records([['bare', 'ACDE'], ['blankdesc', '  ', 'AC-DE']])
    cases = OrderedDict()
    for name, records in (('fasta_parsed', [fasta_rec]),
                          ('clustal_parsed', [clustal_rec]),
                          ('constructed', built)):
        def call(records=records):
            fd, path = tempfile.mkstemp()
            os.close(fd)
            try:
                _share.write_fasta(records, path)
                return read_bytes(path).decode('utf-8')
            finally:
                os.remove(path)
        cases[name] = OrderedDict([('records', records_table(records)),
                                   ('result', outcome(call))])
    result['_share.write_fasta'] = cases
    return result


# ---------------------------------------------------------------------------
# Library: update_fasta

UNGAP_CASES = OrderedDict([
    # Known wrong: the X-containing record ends the loop, so the two unique
    # shorter records are lost.
    ('x_redundant_then_unique', [
        ['long', 'ACDEFGHIKL'], ['withX_contained', 'CDE-XGH'],
        ['short_novel1', 'MNP'], ['short_novel2', 'QRS']]),
    ('x_novel_segment', [
        ['long', 'ACDEFGHIKL'], ['withX_novel', 'CDEXWWW'],
        ['short_novel', 'MNP']]),
    ('x_segments_in_different_records', [
        ['first', 'AAAAACCCCC'], ['second', 'DDDDDEEEEE'],
        ['split', 'CCCXDDD'], ['tail', 'MNP']]),
    ('exact_duplicates', [
        ['a', 'ACDEF'], ['b', 'ACDEF'], ['c', 'GHIK']]),
    ('contained_substring', [
        ['full', 'MACDEFGHIK'], ['frag', 'DEFG'], ['other', 'WWYY']]),
    ('gaps_dots_and_case', [
        ['gappy', 'ac-de.fg'], ['plain', 'ACDEFG'], ['dots', '..KL..']]),
    ('length_ties', [
        ['b_id', 'MNPQ'], ['a_id', 'ACDE'], ['c_id', 'WXYZ']]),
    ('empty_and_all_gap', [
        ['allgap', '----'], ['real', 'ACDE'], ['empty', '']]),
])


def capture_update_fasta():
    from Bio import SeqIO
    from fammerlib import update_fasta
    result = OrderedDict()
    cases = OrderedDict()
    for name, rows in UNGAP_CASES.items():
        cases[name] = OrderedDict([
            ('input', rows),
            ('result', with_log(lambda rows=rows: records_table(
                update_fasta.ungap_and_unique(make_records(rows),
                                              'target.fasta')))),
        ])
    pkl = os.path.join(REPO, 'profiles', 'ROPK', 'PKL.aln')
    cases['profiles/ROPK/PKL.aln'] = OrderedDict([
        ('input', 'profiles/ROPK/PKL.aln'),
        ('result', with_log(lambda: records_table(
            update_fasta.ungap_and_unique(SeqIO.parse(pkl, 'clustal'),
                                          'PKL.fasta')))),
    ])
    result['ungap_and_unique'] = cases

    aln_text = ('CLUSTAL X (-like) multiple sequence alignment\n\n'
                's1 AC-DE\ns2 A--DE\ns3 acdef\n')
    fasta_cases = OrderedDict([
        ('same', '>s1\nACDE\n>s2\nADE\n>s3\nACDEF\n'),
        ('same_other_order', '>s3\nACDEF\n>s1\nACDE\n>s2\nADE\n'),
        ('lowercase', '>s1\nacde\n>s2\nade\n>s3\nacdef\n'),
        ('with_description', '>s1 desc\nACDE\n>s2\nADE\n>s3\nACDEF\n'),
        ('gapped_fasta', '>s1\nAC-DE\n>s2\nADE\n>s3\nACDEF\n'),
        ('extra_sequence', '>s1\nACDE\n>s2\nADE\n>s3\nACDEF\n>s4\nWW\n'),
        ('missing_sequence', '>s1\nACDE\n>s2\nADE\n'),
    ])
    workdir = tempfile.mkdtemp(prefix='fammer-update-')
    try:
        aln_path = os.path.join(workdir, 'x.aln')
        write_bytes(aln_path, aln_text.encode('ascii'))
        same = OrderedDict()
        for name, text in fasta_cases.items():
            fa_path = os.path.join(workdir, 'x.fasta')
            write_bytes(fa_path, text.encode('ascii'))
            same[name] = OrderedDict([
                ('fasta', text),
                ('result', outcome(update_fasta.is_same_aln_fasta,
                                   aln_path, fa_path)),
            ])
        result['is_same_aln_fasta'] = OrderedDict([('aln', aln_text),
                                                   ('cases', same)])
    finally:
        shutil.rmtree(workdir)
    return result


# ---------------------------------------------------------------------------
# Library: tmalign

def trim_cpu_time(text):
    """TM-align output without its final '#Total CPU time' line."""
    lines = text.splitlines(True)
    while lines and lines[-1].startswith('#Total CPU time'):
        lines.pop()
    return ''.join(lines)


def capture_tmalign(tm_outputs):
    from Bio.SeqRecord import SeqRecord
    from fammerlib import tmalign
    result = OrderedDict()

    def parse(text, ref, eqv):
        ref_rec, eqv_rec, tmscore = tmalign.read_tmalign_as_seqrec_pair(
            text, ref, eqv)
        return OrderedDict([('records', records_table([ref_rec, eqv_rec])),
                            ('tmscore', tmscore)])

    parsed = OrderedDict()
    for relpath, (ref, eqv) in tm_outputs.items():
        text = read_bytes(os.path.join(TOOL_OUT, relpath))
        parsed[relpath] = OrderedDict([
            ('ref_id', ref), ('eqv_id', eqv),
            ('as_written', outcome(parse, text, ref, eqv)),
            # The same output with the final CPU-time line removed, which is
            # the layout the parser was written for.
            ('without_cpu_time_line', outcome(parse, trim_cpu_time(text),
                                              ref, eqv)),
        ])
    result['read_tmalign_as_seqrec_pair'] = parsed

    descriptions = ['TMalign TM-score=0.854630', 'x TM-score=0.5 y',
                    'TM-score=abc', 'no score here', '',
                    'TM-score0.7']
    result['tmscore_from_description'] = [
        [text, outcome(tmalign.tmscore_from_description, text)]
        for text in descriptions]

    dedupe_rows = [
        ['_seed_1ATP.pdb', 'TMalign TM-score=0.800000', 'AC-DE'],
        ['1ATP.pdb', 'TMalign TM-score=0.900000', 'AC-DE'],
        ['_seed_1ATP.pdb', 'TMalign TM-score=0.700000', 'ACD-E'],
        ['1CSN.pdb', 'TMalign TM-score=0.600000', 'WW-YY'],
        ['_seed_1CSN.pdb', 'TMalign TM-score=0.650000', 'WW-YY'],
        ['other', 'TMalign TM-score=0.500000', 'AC-DE'],
    ]
    result['clean_and_dedupe_seqs'] = OrderedDict([
        ('input', dedupe_rows),
        ('best_score_false', outcome(lambda: [
            [rec.id, rec.name, rec.description, str(rec.seq)]
            for rec in tmalign.clean_and_dedupe_seqs(
                make_records(dedupe_rows), best_score=False)])),
        ('best_score_true', outcome(lambda: [
            [rec.id, rec.name, rec.description, str(rec.seq)]
            for rec in tmalign.clean_and_dedupe_seqs(
                make_records(dedupe_rows), best_score=True)])),
    ])
    purge_rows = [['a', 'AC-DE'], ['b', 'AC-DE'], ['c', 'ACD-E'],
                  ['d', 'ac-de']]
    result['purge_seqs'] = OrderedDict([
        ('input', purge_rows),
        ('result', outcome(lambda: records_table(
            tmalign.purge_seqs(make_records(purge_rows))))),
    ])

    # mst_pairs as the code calls it: a fresh SeqRecord pair per TM-align run.
    names = ['A.pdb', 'B.pdb', 'C.pdb', 'D.pdb']
    scores = [0.9, 0.8, 0.7, 0.6, 0.5, 0.4]

    def make_pairs(cls):
        from Bio.Seq import Seq
        pairs = []
        index = 0
        for i, left in enumerate(names):
            for right in names[i + 1:]:
                pairs.append((cls(Seq('ACDE'), id=left),
                              cls(Seq('ACDE'), id=right), scores[index]))
                index += 1
        return pairs

    class IdentityHashedRecord(SeqRecord):
        """SeqRecord hashed by identity, as before Biopython 1.67."""
        __eq__ = object.__eq__
        __ne__ = object.__ne__
        __hash__ = object.__hash__

    def mst(cls):
        return sorted(sorted([left.id, right.id])
                      for left, right in tmalign.mst_pairs(make_pairs(cls)))

    result['mst_pairs'] = OrderedDict([
        ('structures', names),
        ('scores_in_pair_order', scores),
        ('seqrecord', outcome(mst, SeqRecord)),
        # Not the behavior of this environment: an emulation of Biopython
        # before 1.67, recorded for the tmalign seeding decision.
        ('identity_hash_emulation', outcome(mst, IdentityHashedRecord)),
    ])
    return result


# ---------------------------------------------------------------------------
# Library: cluster

def capture_cluster_library():
    from fammerlib import cluster
    result = OrderedDict()
    vectors = OrderedDict([
        ('three_values', [0.1, 0.2, 0.4]),
        ('single', [0.3]),
        ('n100_of_0.05', [0.05] * 100),
        ('n250_of_0.05', [0.05] * 250),
        ('n500_of_0.05', [0.05] * 500),
        ('with_zero', [0.05] * 10 + [0.0]),
        ('empty', []),
    ])
    result['geom_mean'] = OrderedDict(
        (name, OrderedDict([('input', vec if len(vec) < 20 else
                             '%d x %r' % (len(vec), vec[0])),
                            ('result', outcome(cluster.geom_mean, vec))]))
        for name, vec in vectors.items())

    cases = OrderedDict()
    for name, rows, do_ungap in (
            ('with_description', [['s1', 's1 desc', 'AC-de.']], False),
            ('with_description_ungap', [['s1', 's1 desc', 'AC-de.']], True),
            ('empty_description', [['s2', '', 'AC-DE']], False)):
        def call(rows=rows, do_ungap=do_ungap):
            handle = io.BytesIO()
            for rec in make_records(rows):
                cluster.write_fasta(rec, handle, do_ungap=do_ungap)
            return handle.getvalue().decode('utf-8')
        cases[name] = OrderedDict([('input', rows), ('do_ungap', do_ungap),
                                   ('result', outcome(call))])
    result['write_fasta'] = cases
    # cmd_cluster accepts only .aln and .fasta input.
    result['load_tree_other_suffix'] = outcome(cluster.load_tree,
                                               'sequences.txt')
    return result


# ---------------------------------------------------------------------------
# End to end: command-line runs in a scratch directory

class Workspace(object):
    """A scratch directory holding a copy of profiles/SomeKinases and the
    sequence files of test/."""

    def __init__(self, pythonpath, env_path):
        self.root = tempfile.mkdtemp(prefix='fammer-e2e-')
        copy_tracked('profiles/SomeKinases', self.root)
        for fname in ('twenty.fasta', 'bigset.fasta', 'TgME49.pep'):
            shutil.copy(os.path.join(REPO, 'test', fname), self.root)
        self.env = dict(os.environ, PYTHONPATH=pythonpath, PATH=env_path,
                        PYTHONHASHSEED='0')
        NORMALIZE.add(self.root, '<workdir>')

    def run(self, script, argv):
        """Run fammer.py or tmalign.py; record output and changed files."""
        before = snapshot(self.root)
        proc = subprocess.Popen([sys.executable, os.path.join(REPO, script)]
                                + argv, cwd=self.root, env=self.env,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        after = snapshot(self.root)
        return OrderedDict([
            ('script', script),
            ('argv', argv),
            ('returncode', proc.returncode),
            ('stdout', NORMALIZE(stdout.decode('utf-8'))),
            ('stderr', NORMALIZE(stderr.decode('utf-8'))),
            ('created', sorted(set(after) - set(before))),
            ('modified', sorted(name for name in after
                                if name in before
                                and after[name] != before[name])),
            ('deleted', sorted(set(before) - set(after))),
        ])

    def path(self, *parts):
        return os.path.join(self.root, *parts)

    def save_to(self, other_root):
        shutil.rmtree(other_root)
        shutil.copytree(self.root, other_root)

    def restore_from(self, other_root):
        shutil.rmtree(self.root)
        shutil.copytree(other_root, self.root)

    def close(self):
        shutil.rmtree(self.root)


def normalized_hmm(data):
    """An HMMER3 profile without its DATE line."""
    return b''.join(line for line in data.splitlines(True)
                    if not line.startswith(b'DATE '))


def normalized_tblout(data):
    """An hmmsearch --tblout table without the run-specific trailer lines."""
    return b''.join(line for line in data.splitlines(True)
                    if not line.startswith((b'# Date:', b'# Current dir:')))


# File types written by `fammer build` that are kept in full; the rest are
# recorded by checksum only.
KEPT_SUFFIXES = ('.aln', '.cons.seq', '.families.fa', '.nwk')


def file_manifest(ws, names):
    """[name, size, sha256, kind] rows; .hmm hashed without its DATE line.

    hmmpress's binary files embed the build date, so only their size is
    recorded."""
    rows = []
    for name in sorted(names):
        data = read_bytes(ws.path(name))
        if name.endswith('.hmm'):
            rows.append([name, len(data), sha256(normalized_hmm(data)),
                         'sha256 of the file without its DATE line'])
        elif re.search(r'\.h3[fimp]$', name):
            rows.append([name, len(data), None, 'size only'])
        elif name.endswith('.tbl'):
            normalized = normalized_tblout(data)
            rows.append([name, len(normalized), sha256(normalized),
                         'size and sha256 without the # Date and '
                         '# Current dir lines'])
        else:
            rows.append([name, len(data), sha256(data), 'sha256'])
    return rows


def capture_e2e(ws):
    from fammerlib import _share, add
    e2e = OrderedDict()

    # 1. build --hmmer --no-pdb --tree
    run = ws.run('fammer.py', ['build', 'SomeKinases', '--hmmer', '--no-pdb',
                               '--tree'])
    run['manifest'] = file_manifest(ws, run['created'] + run['modified'])
    kept_root = os.path.join(OUT, 'e2e', 'build')
    for name in run['created']:
        if name.endswith(KEPT_SUFFIXES):
            write_bytes(os.path.join(kept_root, name),
                        read_bytes(ws.path(name)))
    run['kept_files'] = ('e2e/build/: every created file ending in '
                         + ', '.join(KEPT_SUFFIXES))
    e2e['build'] = run
    hmm_files = sorted(name for name in run['created']
                       if name.endswith('.hmm') and '_all' not in name)

    # 2. Functions that call hmmsearch, and HMM-file helpers, on the build.
    cwd = os.getcwd()
    os.chdir(ws.root)
    saved_path = os.environ['PATH']
    os.environ['PATH'] = ws.env['PATH']
    try:
        lib = OrderedDict()
        lib['get_hmm_evalues'] = OrderedDict(
            (name, outcome(lambda name=name: sorted(_share.get_hmm_evalues(
                name, name[:-len('.hmm')] + '.fasta').items())))
            for name in ('SomeKinases/CMGC/CDK.hmm', 'SomeKinases/Aurora.hmm',
                         'SomeKinases/STE/STE7.hmm'))
        lib['get_hmm_len'] = [[name, outcome(add.get_hmm_len, name)]
                              for name in hmm_files]
        lib['walk_profiles'] = outcome(lambda: sorted(
            list(pair) for pair in add.walk_profiles('SomeKinases')))
        lib['note'] = ('Run in the e2e workspace after `fammer build --hmmer '
                       '--no-pdb`; paths are relative to it.')
    finally:
        os.environ['PATH'] = saved_path
        os.chdir(cwd)
    write_json(os.path.join('library', 'hmmer.json'), lib)

    # 3. scan
    scans = OrderedDict()
    for name, argv in (
            ('twenty.summary', ['scan', 'SomeKinases', 'twenty.fasta']),
            ('twenty.table', ['scan', 'SomeKinases', 'twenty.fasta',
                              '--table']),
            ('twenty.include', ['scan', 'SomeKinases', 'twenty.fasta',
                                '--table', '-n', 'CAMK,CMGC']),
            ('twenty.exclude', ['scan', 'SomeKinases', 'twenty.fasta',
                                '--table', '-x', 'CLK,DYRK']),
            ('twenty.seqout', ['scan', 'SomeKinases_all.hmm', 'twenty.fasta',
                               '-o', 'twenty.hits.fasta', '-O']),
            ('twenty.strict', ['scan', 'SomeKinases', 'twenty.fasta',
                               '--table', '-E', '1e-100']),
            ('TgME49.summary', ['scan', 'SomeKinases', 'TgME49.pep']),
            ('TgME49.table', ['scan', 'SomeKinases', 'TgME49.pep',
                              '--table']),
    ):
        run = ws.run('fammer.py', argv)
        written = OrderedDict()
        for fname in run['created'] + run['modified']:
            if fname.endswith('.tbl'):
                continue
            written[fname] = read_bytes(ws.path(fname)).decode('utf-8')
        run['written_files'] = written
        run['manifest'] = file_manifest(ws, run['created'] + run['modified'])
        scans[name] = run
        if name == 'twenty.table':
            tbl = normalized_tblout(read_bytes(ws.path('twenty.fasta.tbl')))
            write_bytes(os.path.join(TOOL_OUT, 'hmmsearch',
                                     'twenty.fasta.SomeKinases_all.tbl'), tbl)
    e2e['scan'] = scans

    # 4. update-fasta, on a copy of the built tree
    built = tempfile.mkdtemp(prefix='fammer-built-')
    ws.save_to(built)
    upd = OrderedDict()
    for name, argv, edited in (
            ('tree', ['update-fasta', 'SomeKinases'], None),
            ('single_aln', ['update-fasta', 'SomeKinases/CMGC/CDK.aln'],
             None),
            ('single_aln_after_edit',
             ['update-fasta', 'SomeKinases/CMGC/CDK.aln'],
             'SomeKinases/CMGC/CDK.aln'),
            ('tree_after_edit', ['update-fasta', 'SomeKinases'],
             'SomeKinases/CMGC/DYRK.aln'),
    ):
        if edited:
            # A curator deletes the last record from an alignment.
            lines = read_bytes(ws.path(edited)).splitlines(True)
            write_bytes(ws.path(edited), b''.join(lines[:-1]))
        run = ws.run('fammer.py', argv)
        if edited:
            run['edit'] = ('Removed the last line, one record, of the built '
                           + edited)
        run['changed_files'] = OrderedDict(
            (fname, read_bytes(ws.path(fname)).decode('utf-8'))
            for fname in run['modified'] + run['created'])
        upd[name] = run
    e2e['update-fasta'] = upd

    # 5. add, starting from the built tree
    ws.restore_from(built)
    sizes = dict((name, len(read_bytes(ws.path(name))))
                 for name in snapshot(ws.root) if name.endswith('.fasta'))
    run = ws.run('fammer.py', ['add', 'SomeKinases', 'TgME49.pep'])
    run['appended'] = OrderedDict(
        (fname, read_bytes(ws.path(fname))[sizes.get(fname, 0):]
         .decode('utf-8'))
        for fname in run['modified'] + run['created']
        if fname.endswith('.fasta'))
    run['manifest'] = file_manifest(ws, run['created'] + run['modified'])
    e2e['add'] = run

    # 6. refine --dry-run, starting from the built tree
    ws.restore_from(built)
    run = ws.run('fammer.py', [
        'refine', 'SomeKinases', '-n',
        '-u', 'SomeKinases/Testfam-Unclassified.fasta'])
    e2e['refine'] = run

    # 6b. build --clean on the built tree: nothing to rebuild, so the run
    #     only deletes each task's intermediate files.
    ws.restore_from(built)
    run = ws.run('fammer.py', ['build', 'SomeKinases', '--hmmer', '--no-pdb',
                               '--clean'])
    e2e['build-clean'] = run

    # 7. The structure-alignment path, which fails on this toolchain.  The
    #    --no-pdb build left empty .pdb.seq files newer than the PDB files;
    #    remove them so that this build aligns the structures.
    ws.restore_from(built)
    for name in snapshot(ws.root):
        if name.endswith('.pdb.seq'):
            os.remove(ws.path(name))
    pdb = OrderedDict()
    run = ws.run('fammer.py', ['build', 'SomeKinases', '--hmmer'])
    run['manifest'] = file_manifest(ws, run['created'] + run['modified'])
    pdb['build_with_pdb'] = run
    ws_cmgc = 'SomeKinases/CMGC/'
    for name, argv in (
            ('tmalign.two', [ws_cmgc + '1JST.pdb', ws_cmgc + '2ERK.pdb']),
            ('tmalign.three', [ws_cmgc + '1JST.pdb', ws_cmgc + '2ERK.pdb',
                               ws_cmgc + '1O9U.pdb']),
    ):
        run = ws.run('tmalign.py', argv + ['-o', 'tm.out.fasta'])
        pdb[name] = run
    e2e['pdb'] = pdb
    shutil.rmtree(built)

    # 8. cluster, run in-process so the MAFFT and FastTree outputs it reads
    #    can be kept as tool outputs.
    e2e['cluster'] = capture_cluster_runs(ws)
    return e2e


def capture_cluster_runs(ws):
    from fammerlib import cluster
    runs = OrderedDict()
    real_check_output = subprocess.check_output
    cwd = os.getcwd()
    saved_path = os.environ['PATH']
    os.environ['PATH'] = ws.env['PATH']
    try:
        for case, seqfile, label in (
                # Identifiers such as 'TAIR|locus:2011201' contain ':'.
                ('CDK.aln', ws.path('SomeKinases/CMGC/CDK.aln'),
                 'SomeKinases/CMGC/CDK.aln, as built above'),
                # Unaligned input: cluster aligns it with mafft --auto.
                ('DYRK.fasta', ws.path('SomeKinases/CMGC/DYRK.fasta'),
                 'profiles/SomeKinases/CMGC/DYRK.fasta'),
                # 631 sequences: the geometric mean of the tip branch
                # lengths underflows.
                ('PKL.aln', os.path.join(REPO, 'profiles', 'ROPK', 'PKL.aln'),
                 'profiles/ROPK/PKL.aln'),
        ):
            workdir = tempfile.mkdtemp(prefix='fammer-cluster-')
            NORMALIZE.add(workdir, '<clusterdir>')
            shutil.copy(seqfile, workdir)
            seqname = os.path.basename(seqfile)
            calls = []

            def recording_check_output(argv, *args, **kwargs):
                data = real_check_output(argv, *args, **kwargs)
                calls.append((list(argv), data))
                return data

            os.chdir(workdir)
            cluster.subprocess.check_output = recording_check_output
            try:
                before = set(os.listdir(workdir))
                namespace = argparse.Namespace(sequences=seqname)
                result = with_log(cluster.cmd_cluster, namespace)
            finally:
                cluster.subprocess.check_output = real_check_output
                os.chdir(cwd)
            tool_files = []
            for index, (argv, data) in enumerate(calls):
                program = os.path.basename(argv[0])
                rel = os.path.join('cluster', '%s.%d.%s.out'
                                   % (case, index, program))
                write_bytes(os.path.join(TOOL_OUT, rel), data)
                tool_files.append([[NORMALIZE(arg) for arg in argv],
                                   'tool-outputs/' + rel])
            # Clusters of equal size are numbered in dictionary order of
            # Bio.Phylo Clade objects, which Python 2 hashes by identity, so
            # the numbering (and the colors in the phyloXML file) can differ
            # between runs.  Record the clusters without their numbers.
            written = sorted(set(os.listdir(workdir)) - before)
            clusters = []
            unique = None
            for fname in written:
                if fname.endswith('.xml'):
                    continue
                headers = [line[1:] for line in read_bytes(
                    os.path.join(workdir, fname)).decode('utf-8').splitlines()
                    if line.startswith('>')]
                if fname.endswith('.Unique'):
                    unique = headers
                else:
                    clusters.append(headers)
            clusters.sort(key=lambda headers: (-len(headers), headers))
            outputs = OrderedDict([
                ('files_written', written),
                ('clusters', clusters),
                ('unique', unique),
            ])
            for row in result['log']:
                row[2] = NORMALIZE(row[2])
            runs[case] = OrderedDict([
                ('input', label),
                ('result', result),
                ('tool_calls', tool_files),
                ('outputs', outputs),
            ])
            shutil.rmtree(workdir)
    finally:
        os.environ['PATH'] = saved_path
    return runs


def capture_tm_outputs(env_path):
    """Run TM-align on fixed pairs; keep its output with the CPU time zeroed."""
    kept = OrderedDict()
    env = dict(os.environ, PATH=env_path)
    for ref, eqv in TMALIGN_PAIRS:
        proc = subprocess.Popen(['TMalign', ref, eqv], cwd=REPO, env=env,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT)
        output = proc.communicate()[0].decode('utf-8')
        assert proc.returncode == 0, output
        output = re.sub(r'#Total CPU time is +[0-9.]+ seconds',
                        '#Total CPU time is  0.00 seconds', output)
        rel = os.path.join('tmalign', '%s__%s.txt' % (
            os.path.basename(ref)[:-4], os.path.basename(eqv)[:-4]))
        write_bytes(os.path.join(TOOL_OUT, rel), output.encode('utf-8'))
        kept[rel] = (ref, eqv)
    return kept


# ---------------------------------------------------------------------------

def main():
    assert sys.version_info[:2] == (2, 7), sys.version
    assert not sys.flags.hash_randomization, "unset PYTHONHASHSEED or -R"
    missing = [tool for tool in TOOLS
               if not any(os.access(os.path.join(d, tool), os.X_OK)
                          for d in os.environ['PATH'].split(os.pathsep))]
    assert not missing, "not on PATH: %s" % ', '.join(missing)
    def last_commit(paths):
        return subprocess.check_output(
            ['git', 'log', '-1', '--format=%H', '--'] + paths,
            cwd=REPO).decode('ascii').strip()
    dirty = subprocess.check_output(['git', 'status', '--porcelain', '--']
                                    + SOURCES + INPUTS,
                                    cwd=REPO).decode('ascii')
    assert not dirty.strip(), "fammer sources or inputs are modified:\n" + dirty
    assert last_commit(SOURCES) == FAMMER_SOURCES_COMMIT, (
        "fammer sources are not those of %s; run from a worktree as "
        "README.md describes" % FAMMER_SOURCES_COMMIT[:7])

    upstream = tempfile.mkdtemp(prefix='fammer-py2-upstream-')
    NORMALIZE.add(upstream, '<upstream>')
    NORMALIZE.add(REPO, '<repo>')
    NORMALIZE.add(sys.prefix, '<python-prefix>')
    NORMALIZE.add(tempfile.gettempdir(), '<tmpdir>')
    siblings = os.path.dirname(REPO)
    git_archive(os.environ.get('BIOFRILLS_REPO',
                               os.path.join(siblings, 'biofrills')),
                BIOFRILLS_COMMIT, 'biofrills', upstream)
    git_archive(os.environ.get('BIOCMA_REPO',
                               os.path.join(siblings, 'biocma')),
                BIOCMA_COMMIT, 'biocma', upstream)
    sys.path[:0] = [REPO, upstream]
    pythonpath = REPO + os.pathsep + upstream

    import Bio
    import networkx
    import reportlab
    import biofrills.consensus
    import biocma.biocma
    assert biofrills.consensus.__file__.startswith(upstream)
    assert biocma.biocma.__file__.startswith(upstream)

    for path in (OUT, TOOL_OUT):
        if os.path.isdir(path):
            shutil.rmtree(path)
    try:
        env_path = SHIMS + os.pathsep + os.environ['PATH']
        tm_outputs = capture_tm_outputs(env_path)
        ws = Workspace(pythonpath, env_path)
        try:
            e2e = capture_e2e(ws)
        finally:
            ws.close()
        for name, data in e2e.items():
            write_json(os.path.join('e2e', name + '.json'), data)

        write_json(os.path.join('library', 'tasks.json'), capture_tasks())
        for tree in ('SomeKinases', 'ROPK'):
            write_json(os.path.join('library', 'build-%s.json' % tree),
                       capture_build_graph(tree))
        write_json(os.path.join('library', 'build-SomeKinases-group-fasta.json'),
                   capture_build_graph('SomeKinases', group_fastas=True))
        write_json(os.path.join('library', 'mapgaps.json'),
                   capture_mapgaps(os.path.join(OUT, 'e2e', 'build')))
        write_json(os.path.join('library', 'share-scan-add.json'),
                   capture_share(os.path.join(
                       TOOL_OUT, 'hmmsearch',
                       'twenty.fasta.SomeKinases_all.tbl')))
        write_json(os.path.join('library', 'update_fasta.json'),
                   capture_update_fasta())
        write_json(os.path.join('library', 'tmalign.json'),
                   capture_tmalign(tm_outputs))
        write_json(os.path.join('library', 'cluster.json'),
                   capture_cluster_library())
    finally:
        shutil.rmtree(upstream)

    write_json('versions.json', OrderedDict([
        ('python', sys.version.split()[0]),
        ('implementation', platform.python_implementation()),
        ('platform', platform.platform()),
        ('machine', platform.machine()),
        ('biopython', Bio.__version__),
        ('networkx', networkx.__version__),
        ('reportlab', reportlab.Version),
        ('biofrills_commit', BIOFRILLS_COMMIT),
        ('biocma_commit', BIOCMA_COMMIT),
        ('fammer_sources_commit', FAMMER_SOURCES_COMMIT),
        ('fammer_inputs_commit', last_commit(INPUTS)),
        ('mafft', tool_version(['mafft', '--version'], r'v\d')),
        ('hmmer', tool_version(['hmmbuild', '-h'], r'HMMER \d')),
        ('fasttree', tool_version(['fasttree'], r'FastTree.*\d')),
        ('tmalign', tool_version(['TMalign'], r'Version')),
        ('mafft_shim', 'tool-shims/mafft replaces --globalgenafpair with '
                       '--globalpair'),
    ]))


if __name__ == '__main__':
    main()
