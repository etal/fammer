"""The cases recorded under Python 2.7, as calls on the current code.

tests/data/README.md describes the fixtures under tests/data/py27 and
documents every reason a current outcome may differ from a recorded one.
tests/test_reference.py and tests/data/replay.py both use the cases defined
here, so that the test and the differences.json that replay.py writes compare
outcomes in the same way.

A library case is a call recorded in one of the files under py27/library.
Its outcome is compared with the recorded one exactly, with integers and
floats kept distinct, except that an exception is compared by its type alone,
since the message text is the interpreter's, and a recorded built-in
exception type also accepts its subclasses, so that the FileNotFoundError of
Python 3 matches the OSError that Python 2 raised.  Where the recorded outcome
lists the messages logged at level INFO and above, or the text printed to
standard output, the current call is recorded in the same way.

The update-fasta runs of py27/e2e/update-fasta.json are replayed in-process,
in order, on a scratch copy of profiles/SomeKinases with the alignments that
the recorded `fammer build` wrote, and are compared as described in
update_fasta_run_differences.
"""

import builtins
import contextlib
import hashlib
import importlib.util
import io
import json
import logging
import os
import re
import shutil
import subprocess
import tempfile
from collections.abc import Callable
from dataclasses import dataclass
from functools import cache, partial
from pathlib import Path

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biofrills import consensus

from fammer import _share, add, build, cluster, scan, tasks, tmalign, update_fasta

TESTS = Path(__file__).resolve().parent
REPO = TESTS.parent
DATA = TESTS / "data"
FIXTURES = DATA / "py27"
LIBRARY = FIXTURES / "library"
TOOL_OUTPUTS = DATA / "tool-outputs"

# Requirements a case may name, besides the current code: the input
# directories profiles/ and test/, which an sdist does not contain and which
# are copied as Git tracks them, and the optional dependency biocma.
CHECKOUT = "checkout"
BIOCMA = "biocma"


@dataclass
class Case:
    """One recorded call: where it is recorded, how to repeat it, its outcome.

    ``case`` is the name or the position of the call in its section.  ``call``
    returns the value as the fixture records it.  ``requires`` names what the
    call needs besides the current code: CHECKOUT, BIOCMA or both.
    """

    file: str
    section: str
    case: str | int
    call: Callable
    recorded: dict
    requires: tuple = ()

    @property
    def key(self):
        return (self.file, self.section, self.case)


def entry_key(entry):
    """The Case.key of the case that a differences.json entry names."""
    return (entry["file"], entry["section"], entry["case"])


def is_checkout():
    """Whether the repository is a Git checkout with its input directories."""
    return (REPO / ".git").exists()


def missing_requirements(case):
    """The requirements of a case that this environment does not meet."""
    missing = []
    for name in case.requires:
        if name == CHECKOUT and not is_checkout():
            missing.append(name)
        elif name == BIOCMA and importlib.util.find_spec("biocma") is None:
            missing.append(name)
    return missing


@cache
def load(relpath):
    """A fixture file, by its path relative to py27/."""
    with open(FIXTURES / relpath, encoding="ascii") as handle:
        return json.load(handle)


# ---------------------------------------------------------------------------
# Outcomes and their comparison


class LogRecorder(logging.Handler):
    """Collect log records as [level, module, message] rows."""

    def __init__(self):
        super().__init__()
        self.rows = []

    def emit(self, record):
        self.rows.append([record.levelname, record.module, record.getMessage()])


def outcome(func, *args, log=False, stdout=False, **kwargs):
    """Call func and return its outcome as the capture records it.

    The value goes through a JSON round trip, which turns tuples into lists
    and reads floats back exactly.  With ``log``, the messages logged at
    level INFO and above are added as ``log``; with ``stdout``, the printed
    text as ``stdout``.
    """
    recorder = LogRecorder()
    root = logging.getLogger()
    saved_level = root.level
    printed = io.StringIO()
    if log:
        root.setLevel(logging.INFO)
        root.addHandler(recorder)
    try:
        with (
            contextlib.redirect_stdout(printed) if stdout else contextlib.nullcontext()
        ):
            try:
                value = func(*args, **kwargs)
            except Exception as exc:
                result = {"error": {"type": type(exc).__name__}}
            else:
                result = {"value": json.loads(json.dumps(value))}
    finally:
        if log:
            root.removeHandler(recorder)
            root.setLevel(saved_level)
    if log:
        result["log"] = recorder.rows
    if stdout:
        result["stdout"] = printed.getvalue()
    return result


def case_outcome(case):
    """The current outcome of a case, recorded as its fixture records it."""
    return outcome(
        case.call, log="log" in case.recorded, stdout="stdout" in case.recorded
    )


def comparable(recorded):
    """A recorded outcome without the message text of any exception in it."""
    if isinstance(recorded, dict):
        return {
            key: (
                {"type": value["type"]}
                if key == "error" and isinstance(value, dict)
                else comparable(value)
            )
            for key, value in recorded.items()
        }
    if isinstance(recorded, list):
        return [comparable(value) for value in recorded]
    return recorded


def same_error_type(actual, expected):
    """Whether an exception type name matches a recorded one.

    A recorded built-in type also accepts the built-in types derived from it.
    """
    if actual == expected:
        return True
    actual_cls = getattr(builtins, actual, None)
    expected_cls = getattr(builtins, expected, None)
    return (
        isinstance(actual_cls, type)
        and isinstance(expected_cls, type)
        and issubclass(actual_cls, BaseException)
        and issubclass(actual_cls, expected_cls)
    )


def same(actual, expected):
    """Whether a current outcome, or a part of one, agrees with a recorded one."""
    if type(actual) is not type(expected):
        return False
    if isinstance(expected, dict):
        if actual.keys() != expected.keys():
            return False
        for key, value in expected.items():
            if key == "error" and isinstance(value, dict):
                if not same_error_type(actual[key].get("type"), value["type"]):
                    return False
            elif not same(actual[key], value):
                return False
        return True
    if isinstance(expected, list):
        return len(actual) == len(expected) and all(
            same(a, e) for a, e in zip(actual, expected, strict=True)
        )
    return actual == expected


def matches(actual, recorded):
    """Whether a current outcome agrees with a recorded one."""
    return same(actual, comparable(recorded))


# ---------------------------------------------------------------------------
# Helpers shared by the calls


def records_table(records):
    """SeqRecords as [id, description, sequence] rows, as captured."""
    return [[rec.id, rec.description, str(rec.seq)] for rec in records]


def make_records(rows):
    """[id, sequence] or [id, description, sequence] rows as SeqRecords."""
    records = []
    for row in rows:
        if len(row) == 2:
            records.append(SeqRecord(Seq(row[1]), id=row[0]))
        else:
            records.append(SeqRecord(Seq(row[2]), id=row[0], description=row[1]))
    return records


def sha256(data):
    return hashlib.sha256(data).hexdigest()


@cache
def tracked_files(relpath):
    """The Git-tracked files under a repository-relative directory."""
    listing = subprocess.run(
        ["git", "ls-files", "-z", "--", relpath],
        cwd=REPO,
        stdout=subprocess.PIPE,
        text=True,
        check=True,
    )
    names = tuple(name for name in listing.stdout.split("\0") if name)
    if not names:
        raise RuntimeError(f"git lists no files under {relpath}")
    return names


def copy_tracked(relpath, destdir):
    """Copy the Git-tracked files under relpath into destdir, keeping their
    paths relative to the parent directory of relpath."""
    parent = os.path.dirname(relpath)
    for name in tracked_files(relpath):
        dest = os.path.join(destdir, os.path.relpath(name, parent))
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        shutil.copy(REPO / name, dest)


# ---------------------------------------------------------------------------
# library/tasks.json


def run_task_scenario(name):
    """One Task.build scenario of capture_task_scenarios, in a scratch
    directory, as the capture records it."""
    t_old, t_mid, t_new = 1000000000, 1000000100, 1000000200
    calls = []

    with tempfile.TemporaryDirectory(prefix="fammer-tasks-") as workdir:

        def path(fname):
            return os.path.join(workdir, fname)

        def make(fname, mtime):
            with open(path(fname), "wb") as handle:
                handle.write(b"content\n")
            os.utime(path(fname), (mtime, mtime))

        def action(task):
            calls.append(os.path.basename(task.target))
            with open(task.target, "wb") as handle:
                handle.write(b"built\n")

        def failing_action(task):
            calls.append(os.path.basename(task.target))
            raise ValueError("action failed")

        def relative(result):
            for row in result.get("log", []):
                row[2] = row[2].replace(workdir + os.sep, "")
            return result

        if name == "clean":
            for fname in ("a.tmp", "b.tmp", "c.tmp", "out", "mid"):
                make(fname, t_old)
            top = tasks.Task(
                path("out"),
                cleans=[path("a.tmp"), path("absent.tmp")],
                depends=[
                    tasks.Task(path("mid"), cleans=path("b.tmp")),
                    path("c.tmp"),
                ],
            )
            return {
                "clean": relative(outcome(top.clean, log=True)),
                "remaining_files": sorted(os.listdir(workdir)),
            }

        if name == "target_missing":
            make("dep.in", t_old)
            top = tasks.Task(path("out"), action, depends=path("dep.in"))
        elif name == "target_newer_than_dep":
            make("dep.in", t_old)
            make("out", t_new)
            top = tasks.Task(path("out"), action, depends=path("dep.in"))
        elif name == "target_same_mtime_as_dep":
            make("dep.in", t_mid)
            make("out", t_mid)
            top = tasks.Task(path("out"), action, depends=path("dep.in"))
        elif name == "target_older_than_dep":
            make("dep.in", t_new)
            make("out", t_old)
            top = tasks.Task(path("out"), action, depends=path("dep.in"))
        elif name == "stale_middle_task":
            make("dep.in", t_mid)
            make("mid", t_old)
            make("top", t_new)
            mid = tasks.Task(path("mid"), action, depends=path("dep.in"))
            top = tasks.Task(path("top"), action, depends=mid)
        elif name == "missing_file_dependency_with_target":
            make("out", t_new)
            top = tasks.Task(path("out"), action, depends=path("absent.in"))
        elif name == "missing_file_dependency_without_target":
            top = tasks.Task(path("out"), action, depends=path("absent.in"))
        elif name == "failing_action":
            make("dep.in", t_old)
            top = tasks.Task(path("out"), failing_action, depends=path("dep.in"))
        else:
            raise KeyError(name)
        up_to_date = outcome(top.is_up_to_date)
        built = relative(outcome(top.build, log=True))
        return {
            "is_up_to_date_before": up_to_date,
            "build": built,
            "actions_called": list(calls),
            "is_up_to_date_after": outcome(top.is_up_to_date),
        }


def is_empty_of_size(size):
    """tasks.is_empty on a scratch file of the given size in bytes."""
    with tempfile.TemporaryDirectory(prefix="fammer-tasks-") as workdir:
        path = os.path.join(workdir, f"size{size}")
        with open(path, "wb") as handle:
            handle.write(b"x" * size)
        return tasks.is_empty(path)


def tasks_cases():
    file = "library/tasks.json"
    data = load(file)
    for section in ("ext", "noext"):
        for index, (path, recorded) in enumerate(data[section]):
            args = (path, "aln") if section == "ext" else (path,)
            call = partial(getattr(tasks, section), *args)
            yield Case(file, section, index, call, recorded)
    task = tasks.Task("group/fam.aln")
    yield Case(
        file,
        "ext_of_task",
        "group/fam.aln",
        partial(tasks.ext, task, "hmm"),
        data["ext_of_task"],
    )
    yield Case(
        file,
        "noext_of_task",
        "group/fam.aln",
        partial(tasks.noext, task),
        data["noext_of_task"],
    )
    for size, recorded in data["is_empty_by_size"].items():
        call = partial(is_empty_of_size, int(size))
        yield Case(file, "is_empty_by_size", size, call, recorded)
    for name, scenario in data["build_scenarios"].items():
        call = partial(run_task_scenario, name)
        yield Case(file, "build_scenarios", name, call, {"value": scenario})
    for index, (cmd, recorded) in enumerate(data["sh"]):
        yield Case(file, "sh", index, partial(tasks.sh, cmd), recorded)
    call = partial(tasks.which, ["fammer-no-such-program"])
    yield Case(
        file, "which_missing", "fammer-no-such-program", call, data["which_missing"]
    )


# ---------------------------------------------------------------------------
# library/build-*.json


# The group FASTA files that build-SomeKinases-group-fasta.json adds beside
# the group directories, as capture_build_graph writes them.
GROUP_FASTA_CONTENT = b">placeholder\nACDEFGHIK\n"


def describe_dep(dep):
    if isinstance(dep, tasks.Task):
        return {"task": dep.target}
    return {"file": dep}


def task_graph(result):
    """Every task reachable from a build.Result, keyed by target."""
    found = {}
    conflicts = []

    def visit(task):
        if not isinstance(task, tasks.Task):
            return
        entry = {
            "action": task.action.__name__,
            "kwargs": dict(sorted(task.kwargs.items())),
            "depends": [describe_dep(dep) for dep in task.depends],
            "cleans": list(task.cleans),
        }
        if task.target in found:
            if found[task.target] != entry:
                conflicts.append(task.target)
            return
        found[task.target] = entry
        for dep in task.depends:
            visit(dep)

    roots = {}
    for attr in ("aln", "pdbseq", "hmm", "cma", "tpl"):
        task = getattr(result, attr)
        roots[attr] = task.target if isinstance(task, tasks.Task) else task
        visit(task)
    return {
        "result": roots,
        "tasks": dict(sorted(found.items())),
        "conflicting_targets": sorted(set(conflicts)),
    }


def in_tree_copy(relpath, added_files, func, *args):
    """Call func in a scratch directory holding a copy of a profile tree.

    ``relpath`` is the tree's path in the repository, such as
    profiles/SomeKinases, and ``added_files`` are paths of the files that the
    fixture adds to the copy, relative to the scratch directory.
    """
    with tempfile.TemporaryDirectory(prefix="fammer-graph-") as workdir:
        copy_tracked(relpath, workdir)
        for name in added_files:
            with open(os.path.join(workdir, name), "wb") as handle:
                handle.write(GROUP_FASTA_CONTENT)
        with contextlib.chdir(workdir):
            return func(*args)


def sorted_pdb_files(graph):
    """A recorded task graph with the PDB files of each align_pdbs task sorted.

    The Python 2 code listed those files in the order of an unsorted glob,
    which is the enumeration order of the file system, so the recorded order
    is that of the capture's file system; taskify_subdirs now sorts them.  The
    files precede the task dependencies in the list, and only their order is
    normalized.
    """
    graph_tasks = {}
    for target, entry in graph["tasks"].items():
        if entry["action"] == "align_pdbs":
            files = sorted(
                (dep for dep in entry["depends"] if "file" in dep),
                key=lambda dep: dep["file"],
            )
            others = [dep for dep in entry["depends"] if "file" not in dep]
            entry = {**entry, "depends": files + others}
        graph_tasks[target] = entry
    return {**graph, "tasks": graph_tasks}


def taskify_graph(tree, kwargs):
    return task_graph(
        build.taskify_subdirs(
            tree,
            kwargs["hmmer"],
            kwargs["mapgaps"],
            kwargs["use_pdb"],
            kwargs["level"],
        )
    )


BUILD_FIXTURES = [
    "library/build-SomeKinases.json",
    "library/build-ROPK.json",
    "library/build-SomeKinases-group-fasta.json",
]


def build_cases(file):
    data = load(file)
    relpath = data["tree"]
    tree = os.path.basename(relpath)
    added = data.get("added_files", [])
    call = partial(in_tree_copy, relpath, added, build.treeify_subdirs, tree)
    yield Case(
        file, "treeify_subdirs", tree, call, data["treeify_subdirs"], (CHECKOUT,)
    )
    for index, graph in enumerate(data["taskify_subdirs"]):
        call = partial(
            in_tree_copy, relpath, added, taskify_graph, tree, graph["kwargs"]
        )
        recorded = graph["graph"]
        if "value" in recorded:
            recorded = {**recorded, "value": sorted_pdb_files(recorded["value"])}
        yield Case(file, "taskify_subdirs", index, call, recorded, (CHECKOUT,))


# ---------------------------------------------------------------------------
# library/mapgaps.json


def mapgaps_name(target):
    """The profile name that mg_aln2cma gives the consensus of a target."""
    return os.path.basename(tasks.noext(target))


def written_profile(target, write):
    """Call write(path) for a scratch path named like target; return the
    written file's SHA-256 and the paths created under the scratch directory."""
    with tempfile.TemporaryDirectory(prefix="fammer-mapgaps-") as workdir:
        path = os.path.join(workdir, target)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        write(path)
        created = sorted(
            os.path.relpath(os.path.join(dirpath, fname), workdir)
            for dirpath, _dirs, fnames in os.walk(workdir)
            for fname in fnames
        )
        with open(path, "rb") as handle:
            return {"sha256": sha256(handle.read()), "created": created}


def run_mg_aln2cma(rec):
    """build.mg_aln2cma as taskify_subdirs defines the task."""

    def write(path):
        task = tasks.Task(
            path,
            action=build.mg_aln2cma,
            depends=str(FIXTURES / rec["alignment"]),
            kwargs={"level": rec["level"]},
        )
        build.mg_aln2cma(task, **task.kwargs)

    return written_profile(rec["target"], write)


def run_biocma_write(rec):
    """The biocma calls of mg_aln2cma, given the recorded consensus.

    This compares biocma's writer with the recorded profile apart from the
    consensus that biofrills computes.
    """
    from biocma import biocma

    def write(path):
        aln = AlignIO.read(FIXTURES / rec["alignment"], "clustal")
        name = mapgaps_name(rec["target"])
        cons_rec = SeqRecord(
            Seq(rec["consensus"]), id=name, description=name + " consensus"
        )
        aln = MultipleSeqAlignment([cons_rec, *aln])
        cmaln = biocma.ChainMultiAlignment(aln, level=rec["level"])
        biocma.write([cmaln], path, do_iron=True)

    return written_profile(rec["target"], write)


def mapgaps_consensus(rec):
    aln = AlignIO.read(FIXTURES / rec["alignment"], "clustal")
    return consensus.consensus(aln, trim_ends=False, gap_threshold=0.6)


def mapgaps_cases():
    file = "library/mapgaps.json"
    for rec in load(file)["cases"]:
        target = rec["target"]
        yield Case(
            file,
            "consensus",
            target,
            partial(mapgaps_consensus, rec),
            {"value": rec["consensus"]},
        )
        with open(LIBRARY / "mapgaps" / target, "rb") as handle:
            kept = {"sha256": sha256(handle.read()), "created": rec["created"]}
        if "value" in rec["mg_aln2cma"]:
            recorded = {"value": kept}
        else:
            recorded = rec["mg_aln2cma"]
        yield Case(
            file,
            "mg_aln2cma",
            target,
            partial(run_mg_aln2cma, rec),
            recorded,
            (BIOCMA,),
        )
        yield Case(
            file,
            "biocma.write",
            target,
            partial(run_biocma_write, rec),
            {"value": kept},
            (BIOCMA,),
        )


# ---------------------------------------------------------------------------
# library/share-scan-add.json

SHARE_TABLE = "tool-outputs/hmmsearch/twenty.fasta.SomeKinases_all.tbl"


def hits_table(hits):
    """{seq: (profile, evalue)} as sorted [seq, profile, evalue] rows."""
    return [[seq, prof, evalue] for seq, (prof, evalue) in sorted(hits.items())]


def table_text(name):
    """The text of an hmmsearch table of share-scan-add.json, by its name."""
    if name.startswith("tool-outputs/"):
        return (DATA / name).read_text(encoding="utf-8")
    return load("library/share-scan-add.json")["inputs"][name]


def scanned_hits(name):
    return _share.parse_scanned_hits(table_text(name).splitlines(True))


def write_fasta_text(records):
    """_share.write_fasta on records, as the text of the written file."""
    with tempfile.TemporaryDirectory(prefix="fammer-share-") as workdir:
        path = os.path.join(workdir, "out.fasta")
        _share.write_fasta(records, path)
        with open(path, encoding="utf-8") as handle:
            return handle.read()


def share_write_fasta_records(name):
    """The records of a _share.write_fasta case, made as the capture made them."""
    if name == "fasta_parsed":
        text = ">seq1 some description\nACDEFGHIKLMNPQRSTVWY\n"
        return [next(SeqIO.parse(io.StringIO(text), "fasta"))]
    if name == "clustal_parsed":
        return [AlignIO.read(REPO / "profiles" / "ROPK" / "PKL.aln", "clustal")[0]]
    if name == "constructed":
        return make_records([["bare", "ACDE"], ["blankdesc", "  ", "AC-DE"]])
    raise KeyError(name)


def share_cases():
    file = "library/share-scan-add.json"
    data = load(file)
    for name, recorded in data["parse_scanned_hits"].items():

        def call(name=name):
            return hits_table(scanned_hits(name))

        yield Case(file, "parse_scanned_hits", name, call, recorded)
    sources = {"synthetic": "synthetic", "real": SHARE_TABLE}
    for section, func in (
        ("scan.write_summary", scan.write_summary),
        ("scan.write_table", scan.write_table),
    ):
        for name, recorded in data[section].items():

            def call(func=func, name=name):
                return func(scanned_hits(sources[name]))

            yield Case(file, section, name, call, recorded)
    for profile, recorded in data["add.filter_profile"].items():

        def call(profile=profile):
            return sorted(add.filter_profile(scanned_hits(SHARE_TABLE), profile))

        yield Case(file, "add.filter_profile", profile, call, recorded)
    profile = data["add.filter_evalue_profile"]
    for index, (threshold, recorded) in enumerate(data["add.filter_evalue"]):

        def call(threshold=threshold):
            hits = sorted(add.filter_profile(scanned_hits(SHARE_TABLE), profile))
            return list(add.filter_evalue(hits, threshold))

        yield Case(file, "add.filter_evalue", index, call, recorded)
    for name, rec in data["_share.write_fasta"].items():

        def call(name=name, rows=rec["records"]):
            records = share_write_fasta_records(name)
            if records_table(records) != rows:
                raise AssertionError(f"the records of {name} differ from the fixture")
            return write_fasta_text(records)

        requires = (CHECKOUT,) if name == "clustal_parsed" else ()
        yield Case(file, "_share.write_fasta", name, call, rec["result"], requires)


# ---------------------------------------------------------------------------
# library/update_fasta.json


def ungap_rows(rows):
    return records_table(
        update_fasta.ungap_and_unique(make_records(rows), "target.fasta")
    )


def ungap_file(relpath):
    target = os.path.basename(relpath)[: -len(".aln")] + ".fasta"
    with open(REPO / relpath, encoding="utf-8") as handle:
        return records_table(
            update_fasta.ungap_and_unique(SeqIO.parse(handle, "clustal"), target)
        )


def same_aln_fasta(aln_text, fasta_text):
    """update_fasta.is_same_aln_fasta on scratch files with the given text."""
    with tempfile.TemporaryDirectory(prefix="fammer-update-") as workdir:
        aln_path = os.path.join(workdir, "x.aln")
        fa_path = os.path.join(workdir, "x.fasta")
        with open(aln_path, "w", encoding="ascii") as handle:
            handle.write(aln_text)
        with open(fa_path, "w", encoding="ascii") as handle:
            handle.write(fasta_text)
        return update_fasta.is_same_aln_fasta(aln_path, fa_path)


def update_fasta_cases():
    file = "library/update_fasta.json"
    data = load(file)
    for name, case in data["ungap_and_unique"].items():
        if isinstance(case["input"], str):
            call = partial(ungap_file, case["input"])
            requires = (CHECKOUT,)
        else:
            call = partial(ungap_rows, case["input"])
            requires = ()
        yield Case(file, "ungap_and_unique", name, call, case["result"], requires)
    aln_text = data["is_same_aln_fasta"]["aln"]
    for name, case in data["is_same_aln_fasta"]["cases"].items():
        call = partial(same_aln_fasta, aln_text, case["fasta"])
        yield Case(file, "is_same_aln_fasta", name, call, case["result"])


# ---------------------------------------------------------------------------
# library/tmalign.json


def trim_cpu_time(text):
    """TM-align output without its final '#Total CPU time' line."""
    lines = text.splitlines(True)
    while lines and lines[-1].startswith("#Total CPU time"):
        lines.pop()
    return "".join(lines)


def parse_tmalign(relpath, variant, ref, eqv):
    text = (TOOL_OUTPUTS / relpath).read_text(encoding="utf-8")
    if variant == "without_cpu_time_line":
        text = trim_cpu_time(text)
    ref_rec, eqv_rec, tmscore = tmalign.read_tmalign_as_seqrec_pair(text, ref, eqv)
    return {"records": records_table([ref_rec, eqv_rec]), "tmscore": tmscore}


def deduped(rows, best_score):
    return [
        [rec.id, rec.name, rec.description, str(rec.seq)]
        for rec in tmalign.clean_and_dedupe_seqs(
            make_records(rows), best_score=best_score
        )
    ]


class IdentityHashedRecord(SeqRecord):
    """SeqRecord hashed by identity, as before Biopython 1.67."""

    __eq__ = object.__eq__
    __ne__ = object.__ne__
    __hash__ = object.__hash__


def mst(names, scores, cls):
    """tmalign.mst_pairs as the code calls it: a fresh record pair per
    TM-align run, scored in pair order."""
    pairs = []
    index = 0
    for i, left in enumerate(names):
        for right in names[i + 1 :]:
            pairs.append(
                (cls(Seq("ACDE"), id=left), cls(Seq("ACDE"), id=right), scores[index])
            )
            index += 1
    return sorted(
        sorted([left.id, right.id]) for left, right in tmalign.mst_pairs(pairs)
    )


def tmalign_cases():
    file = "library/tmalign.json"
    data = load(file)
    for relpath, rec in data["read_tmalign_as_seqrec_pair"].items():
        for variant in ("as_written", "without_cpu_time_line"):
            call = partial(
                parse_tmalign, relpath, variant, rec["ref_id"], rec["eqv_id"]
            )
            yield Case(
                file,
                "read_tmalign_as_seqrec_pair",
                f"{relpath}:{variant}",
                call,
                rec[variant],
            )
    for index, (text, recorded) in enumerate(data["tmscore_from_description"]):
        call = partial(tmalign.tmscore_from_description, text)
        yield Case(file, "tmscore_from_description", index, call, recorded)
    rows = data["clean_and_dedupe_seqs"]["input"]
    for name, best_score in (("best_score_false", False), ("best_score_true", True)):
        call = partial(deduped, rows, best_score)
        recorded = data["clean_and_dedupe_seqs"][name]
        yield Case(file, "clean_and_dedupe_seqs", name, call, recorded)
    rows = data["purge_seqs"]["input"]

    def purge(rows=rows):
        return records_table(tmalign.purge_seqs(make_records(rows)))

    yield Case(file, "purge_seqs", "result", purge, data["purge_seqs"]["result"])
    names = data["mst_pairs"]["structures"]
    scores = data["mst_pairs"]["scores_in_pair_order"]
    for name, cls in (
        ("seqrecord", SeqRecord),
        ("identity_hash_emulation", IdentityHashedRecord),
    ):
        call = partial(mst, names, scores, cls)
        yield Case(file, "mst_pairs", name, call, data["mst_pairs"][name])


# ---------------------------------------------------------------------------
# library/cluster.json

# The geom_mean inputs, as capture_cluster_library defines them; the fixture
# abbreviates the long ones.
GEOM_MEAN_INPUTS = {
    "three_values": [0.1, 0.2, 0.4],
    "single": [0.3],
    "n100_of_0.05": [0.05] * 100,
    "n250_of_0.05": [0.05] * 250,
    "n500_of_0.05": [0.05] * 500,
    "with_zero": [0.05] * 10 + [0.0],
    "empty": [],
}


def cluster_write_fasta(rows, do_ungap):
    handle = io.StringIO()
    for rec in make_records(rows):
        cluster.write_fasta(rec, handle, do_ungap=do_ungap)
    return handle.getvalue()


def cluster_cases():
    file = "library/cluster.json"
    data = load(file)
    for name, rec in data["geom_mean"].items():
        vector = GEOM_MEAN_INPUTS[name]
        shown = vector if len(vector) < 20 else f"{len(vector)} x {vector[0]!r}"
        if shown != rec["input"]:
            raise AssertionError(f"geom_mean input {name} differs from the fixture")
        call = partial(cluster.geom_mean, vector)
        yield Case(file, "geom_mean", name, call, rec["result"])
    for name, rec in data["write_fasta"].items():
        call = partial(cluster_write_fasta, rec["input"], rec["do_ungap"])
        yield Case(file, "write_fasta", name, call, rec["result"])
    call = partial(cluster.load_tree, "sequences.txt")
    yield Case(
        file,
        "load_tree_other_suffix",
        "sequences.txt",
        call,
        data["load_tree_other_suffix"],
    )


# ---------------------------------------------------------------------------
# All library cases

# The top-level keys of each library fixture that hold inputs or notes rather
# than recorded outcomes.
NOT_SECTIONS = {
    "library/share-scan-add.json": {"inputs", "add.filter_evalue_profile"},
    "library/mapgaps.json": {"cases", "kept_files"},
    "library/build-SomeKinases.json": {"tree"},
    "library/build-ROPK.json": {"tree"},
    "library/build-SomeKinases-group-fasta.json": {"tree", "added_files"},
}

# mapgaps.json records its cases in one list; these are the sections that
# mapgaps_cases makes of it.
MAPGAPS_SECTIONS = ["consensus", "mg_aln2cma", "biocma.write"]


def sections(file):
    """The sections of a library fixture, in file order."""
    if file == "library/mapgaps.json":
        return MAPGAPS_SECTIONS
    return [key for key in load(file) if key not in NOT_SECTIONS.get(file, set())]


# The recorded cases of each library fixture, by fixture, in the order they
# are compared.  library/hmmer.json needs HMMER and built profiles;
# tests/test_integration.py compares it.
CASE_SOURCES = {
    "library/tasks.json": tasks_cases,
    **{file: partial(build_cases, file) for file in BUILD_FIXTURES},
    "library/mapgaps.json": mapgaps_cases,
    "library/share-scan-add.json": share_cases,
    "library/update_fasta.json": update_fasta_cases,
    "library/tmalign.json": tmalign_cases,
    "library/cluster.json": cluster_cases,
}
LIBRARY_FIXTURES = list(CASE_SOURCES)


def library_cases(file):
    """The recorded cases of one library fixture, in file order."""
    return CASE_SOURCES[file]()


def all_cases():
    """Every recorded library case."""
    for file in LIBRARY_FIXTURES:
        yield from library_cases(file)


# ---------------------------------------------------------------------------
# e2e/update-fasta.json

UPDATE_FASTA_RUNS = "e2e/update-fasta.json"


def snapshot(topdir):
    """{relative path: sha256} of every file under topdir."""
    result = {}
    for dirpath, _dirs, fnames in os.walk(topdir):
        for fname in fnames:
            path = os.path.join(dirpath, fname)
            with open(path, "rb") as handle:
                result[os.path.relpath(path, topdir)] = sha256(handle.read())
    return result


def fasta_ids(text):
    return [
        line[1:].split(None, 1)[0] for line in text.splitlines() if line.startswith(">")
    ]


def summarize_fasta(recorded, current):
    """The records of two versions of a FASTA file, for an entry."""
    old, new = fasta_ids(recorded), fasta_ids(current)
    return {
        "recorded": {"records": len(old), "sha256": sha256(recorded.encode("utf-8"))},
        "current": {"records": len(new), "sha256": sha256(current.encode("utf-8"))},
        "only_recorded": [i for i in old if i not in new],
        "only_current": [i for i in new if i not in old],
    }


def log_messages(stderr):
    """The lines of fammer's standard error, sorted, without line numbers.

    fammer logs in the format ``module [@line]: message``, and the line
    numbers of the Python 2 sources differ from those of the current ones.
    The lines are sorted because several commands visit files in the
    enumeration order of the file system.
    """
    return sorted(re.sub(r" \[@\d+\]:", ":", line) for line in stderr.splitlines())


def make_update_fasta_workspace(root):
    """Fill root with a copy of profiles/SomeKinases and the recorded
    alignments of the end-to-end build.

    The FASTA files are given an older modification time than the
    alignments, as `fammer build` leaves them.
    """
    copy_tracked("profiles/SomeKinases", root)
    built = FIXTURES / "e2e" / "build"
    for src in built.rglob("*.aln"):
        dest = Path(root) / src.relative_to(built)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(src, dest)
    for dirpath, _dirs, fnames in os.walk(root):
        for fname in fnames:
            mtime = 1000000100 if fname.endswith(".aln") else 1000000000
            os.utime(os.path.join(dirpath, fname), (mtime, mtime))


def update_fasta_run_differences():
    """The differences.json entries of e2e/update-fasta.json.

    The runs are replayed in order, since each starts from the files the
    previous one left.  A run is compared by its outcome, its log lines
    without source line numbers, and the lists of files created, modified and
    deleted; for each file whose contents differ, the entry summarizes the
    records of both versions.  The log lines are compared in sorted order,
    because update-fasta visits the files in the enumeration order of the
    file system.
    """
    import argparse

    data = load(UPDATE_FASTA_RUNS)
    found = []
    with tempfile.TemporaryDirectory(prefix="fammer-replay-e2e-") as root:
        make_update_fasta_workspace(root)
        with contextlib.chdir(root):
            for name, run in data.items():
                if run["argv"][0] != "update-fasta":
                    raise AssertionError(run["argv"])
                if "edit" in run:
                    # A curator deletes the last record from an alignment.
                    edited = run["edit"].rsplit(" ", 1)[-1]
                    with open(edited, "rb") as handle:
                        lines = handle.read().splitlines(True)
                    with open(edited, "wb") as handle:
                        handle.write(b"".join(lines[:-1]))
                before = snapshot(root)
                result = outcome(
                    update_fasta.cmd_update_fasta,
                    argparse.Namespace(basedir=run["argv"][1]),
                    log=True,
                )
                after = snapshot(root)
                current = {
                    "outcome": {key: result[key] for key in result if key != "log"},
                    "log": sorted(
                        f"{module}: {message}"
                        for _level, module, message in result["log"]
                    ),
                    "created": sorted(set(after) - set(before)),
                    "modified": sorted(
                        n for n in after if n in before and after[n] != before[n]
                    ),
                    "deleted": sorted(set(before) - set(after)),
                }
                recorded = {
                    "outcome": (
                        {"value": None}
                        if run["returncode"] == 0
                        else {"error": f"exit status {run['returncode']}"}
                    ),
                    "log": log_messages(run["stderr"]),
                    "created": run["created"],
                    "modified": run["modified"],
                    "deleted": run["deleted"],
                }
                changed = {}
                for path in sorted(
                    set(run["changed_files"])
                    | set(current["created"])
                    | set(current["modified"])
                ):
                    old = run["changed_files"].get(path, "")
                    with open(path, encoding="utf-8") as handle:
                        new = handle.read()
                    if old != new:
                        changed[path] = summarize_fasta(old, new)
                fields = [key for key in recorded if recorded[key] != current[key]]
                if fields or changed:
                    found.append(
                        {
                            "file": UPDATE_FASTA_RUNS,
                            "case": name,
                            "recorded": {key: recorded[key] for key in fields},
                            "current": {key: current[key] for key in fields},
                            "changed_files": changed,
                        }
                    )
    return found
