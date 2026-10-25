#!/usr/bin/env python3
"""Replay recorded reference cases against the current fammer code.

The update-fasta cases recorded under py27/ (library/update_fasta.json and
e2e/update-fasta.json) are run again with the current sources in this
repository, and each case whose outcome differs from the recorded one is
written to differences.json beside this script.  README.md beside it explains
each difference.

Function calls are compared as the capture records them: a value exactly, an
exception by its type only, since the message text is the interpreter's.
The end-to-end cases are replayed in-process on a scratch copy of
profiles/SomeKinases with the alignments that the recorded `fammer build`
wrote (py27/e2e/build/), so no external program is needed; they are compared
by the log lines without source line numbers, the files changed, and the
records of each changed file.

The script requires Python 3 and Biopython.  Usage, from the repository root:

    python tests/data/replay.py
"""

import argparse
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
FIXTURES = os.path.join(HERE, 'py27')
OUT = os.path.join(HERE, 'differences.json')

sys.path.insert(0, REPO)

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from fammer import update_fasta

assert os.path.dirname(os.path.abspath(update_fasta.__file__)) == \
    os.path.join(REPO, 'fammer'), update_fasta.__file__


def load(relpath):
    with open(os.path.join(FIXTURES, relpath), encoding='ascii') as handle:
        return json.load(handle)


def outcome(func, *args, **kwargs):
    """Call func; return {"value": result} or {"error": {"type": ...}}."""
    try:
        return {'value': func(*args, **kwargs)}
    except Exception as exc:
        return {'error': {'type': type(exc).__name__}}


def comparable(result):
    """A recorded outcome, reduced to what the comparison uses."""
    if 'error' in result:
        reduced = {'error': {'type': result['error']['type']}}
    else:
        reduced = {'value': result['value']}
    if 'log' in result:
        reduced['log'] = result['log']
    return reduced


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


def records_table(records):
    """SeqRecords as [id, description, sequence] rows, as captured."""
    return [[rec.id, rec.description, str(rec.seq)] for rec in records]


def make_records(rows):
    """[id, sequence] rows as SeqRecords."""
    return [SeqRecord(Seq(seq), id=seq_id) for seq_id, seq in rows]


def replay_library():
    """Differences in library/update_fasta.json."""
    fname = 'library/update_fasta.json'
    data = load(fname)
    found = []
    for name, case in data['ungap_and_unique'].items():
        if isinstance(case['input'], str):
            path = os.path.join(REPO, case['input'])
            target = os.path.basename(path)[:-len('.aln')] + '.fasta'
            def call(path=path, target=target):
                return records_table(update_fasta.ungap_and_unique(
                    SeqIO.parse(path, 'clustal'), target))
        else:
            def call(rows=case['input']):
                return records_table(update_fasta.ungap_and_unique(
                    make_records(rows), 'target.fasta'))
        current = with_log(call)
        recorded = comparable(case['result'])
        if current != recorded:
            found.append({'file': fname, 'section': 'ungap_and_unique',
                          'case': name, 'recorded': recorded,
                          'current': current})

    workdir = tempfile.mkdtemp(prefix='fammer-replay-')
    try:
        aln_path = os.path.join(workdir, 'x.aln')
        with open(aln_path, 'w', encoding='ascii') as handle:
            handle.write(data['is_same_aln_fasta']['aln'])
        for name, case in data['is_same_aln_fasta']['cases'].items():
            fa_path = os.path.join(workdir, 'x.fasta')
            with open(fa_path, 'w', encoding='ascii') as handle:
                handle.write(case['fasta'])
            current = outcome(update_fasta.is_same_aln_fasta,
                              aln_path, fa_path)
            recorded = comparable(case['result'])
            if current != recorded:
                found.append({'file': fname, 'section': 'is_same_aln_fasta',
                              'case': name, 'recorded': recorded,
                              'current': current})
    finally:
        shutil.rmtree(workdir)
    return found


def snapshot(topdir):
    """{relative path: sha256} of every file under topdir."""
    result = {}
    for dirpath, _dirs, fnames in os.walk(topdir):
        for fname in fnames:
            path = os.path.join(dirpath, fname)
            with open(path, 'rb') as handle:
                digest = hashlib.sha256(handle.read()).hexdigest()
            result[os.path.relpath(path, topdir)] = digest
    return result


def fasta_ids(text):
    return [line[1:].split(None, 1)[0] for line in text.splitlines()
            if line.startswith('>')]


def summarize_fasta(recorded, current):
    """The records of two versions of a FASTA file, for an entry."""
    old, new = fasta_ids(recorded), fasta_ids(current)
    return {
        'recorded': {'records': len(old), 'sha256': hashlib.sha256(
            recorded.encode('utf-8')).hexdigest()},
        'current': {'records': len(new), 'sha256': hashlib.sha256(
            current.encode('utf-8')).hexdigest()},
        'only_recorded': [i for i in old if i not in new],
        'only_current': [i for i in new if i not in old],
    }


def make_workspace():
    """A scratch copy of profiles/SomeKinases with the recorded alignments.

    The FASTA files are given an older modification time than the
    alignments, as `fammer build` leaves them.
    """
    root = tempfile.mkdtemp(prefix='fammer-replay-e2e-')
    names = subprocess.check_output(
        ['git', 'ls-files', '--', 'profiles/SomeKinases'], cwd=REPO,
        text=True).splitlines()
    for name in names:
        dest = os.path.join(root, os.path.relpath(name, 'profiles'))
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        shutil.copy(os.path.join(REPO, name), dest)
    built = os.path.join(FIXTURES, 'e2e', 'build')
    for dirpath, _dirs, fnames in os.walk(built):
        for fname in fnames:
            if fname.endswith('.aln'):
                src = os.path.join(dirpath, fname)
                dest = os.path.join(root, os.path.relpath(src, built))
                os.makedirs(os.path.dirname(dest), exist_ok=True)
                shutil.copy(src, dest)
    for dirpath, _dirs, fnames in os.walk(root):
        for fname in fnames:
            mtime = 1000000100 if fname.endswith('.aln') else 1000000000
            os.utime(os.path.join(dirpath, fname), (mtime, mtime))
    return root


def replay_e2e():
    """Differences in e2e/update-fasta.json, replayed in order."""
    fname = 'e2e/update-fasta.json'
    data = load(fname)
    found = []
    root = make_workspace()
    cwd = os.getcwd()
    try:
        os.chdir(root)
        for name, run in data.items():
            assert run['argv'][0] == 'update-fasta', run['argv']
            if 'edit' in run:
                # A curator deletes the last record from an alignment.
                edited = run['edit'].rsplit(' ', 1)[-1]
                with open(edited, 'rb') as handle:
                    lines = handle.read().splitlines(True)
                with open(edited, 'wb') as handle:
                    handle.write(b''.join(lines[:-1]))
            before = snapshot(root)
            result = with_log(update_fasta.cmd_update_fasta,
                              argparse.Namespace(basedir=run['argv'][1]))
            after = snapshot(root)
            current = {
                'outcome': {key: result[key] for key in result
                            if key != 'log'},
                'log': ['%s: %s' % (module, message)
                        for _level, module, message in result['log']],
                'created': sorted(set(after) - set(before)),
                'modified': sorted(n for n in after if n in before
                                   and after[n] != before[n]),
                'deleted': sorted(set(before) - set(after)),
            }
            recorded = {
                'outcome': ({'value': None} if run['returncode'] == 0
                            else {'error': 'exit status %d'
                                  % run['returncode']}),
                'log': [re.sub(r' \[@\d+\]:', ':', line)
                        for line in run['stderr'].splitlines()],
                'created': run['created'],
                'modified': run['modified'],
                'deleted': run['deleted'],
            }
            changed = {}
            for path in sorted(set(run['changed_files'])
                               | set(current['created'])
                               | set(current['modified'])):
                old = run['changed_files'].get(path, '')
                with open(path, encoding='utf-8') as handle:
                    new = handle.read()
                if old != new:
                    changed[path] = summarize_fasta(old, new)
            fields = [key for key in recorded if recorded[key] != current[key]]
            if fields or changed:
                entry = {'file': fname, 'case': name}
                entry['recorded'] = {key: recorded[key] for key in fields}
                entry['current'] = {key: current[key] for key in fields}
                entry['changed_files'] = changed
                found.append(entry)
    finally:
        os.chdir(cwd)
        shutil.rmtree(root)
    return found


def main():
    """Replay every update-fasta case and rewrite differences.json."""
    found = replay_library() + replay_e2e()
    text = json.dumps({'differences': found}, indent=1,
                      separators=(',', ': '), ensure_ascii=True)
    with open(OUT, 'w', encoding='ascii', newline='\n') as handle:
        handle.write('%s\n' % text)
    print('%d differing cases written to %s' % (len(found), OUT))


if __name__ == '__main__':
    main()
