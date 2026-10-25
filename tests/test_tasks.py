"""Invariants of fammer.tasks.Task over modification-time orderings.

The tests state general properties of ``is_up_to_date``, ``build`` and
``clean`` and check them against a small model: a target is stale exactly when
it is missing or strictly older than some dependency. A target whose
modification time equals a dependency's counts as up to date. The few scenarios
recorded from the Python 2 code base are compared separately, in the reference
tests.
"""

import itertools
import os
from pathlib import Path

import pytest

from fammer.tasks import Task

# Fixed modification times, far in the past, so that a file written by an
# action (which receives the current time) is newer than every fixed one.
BASE = 1_000_000_000

# A dependency's modification time relative to the target's: older, equal, newer.
OFFSETS = (-10, 0, 10)


def make_file(path, mtime=None):
    """Create a small file, optionally with a fixed modification time."""
    Path(path).write_text("content")
    if mtime is not None:
        os.utime(path, (mtime, mtime))
    return str(path)


class Recorder:
    """An action that records which target it built and then writes it."""

    def __init__(self):
        self.calls = []

    def action(self, task):
        self.calls.append(os.path.basename(task.target))
        Path(task.target).write_text("built")

    def task(self, target, depends=None, cleans=None):
        return Task(str(target), action=self.action, depends=depends, cleans=cleans)


class Boom(Exception):
    """Raised by actions that are meant to fail."""


def fail(task):
    raise Boom(task.target)


def is_stale(target_mtime, dep_mtimes):
    """The model: missing, or strictly older than some dependency."""
    return target_mtime is None or any(dep > target_mtime for dep in dep_mtimes)


@pytest.mark.parametrize(
    "target_present", [False, True], ids=["target-missing", "target-present"]
)
@pytest.mark.parametrize(
    "offsets",
    [combo for n in range(3) for combo in itertools.product(OFFSETS, repeat=n)],
    ids=lambda combo: "deps[" + ",".join(f"{o:+d}" for o in combo) + "]",
)
def test_file_dependencies_follow_the_staleness_model(
    tmp_path, target_present, offsets
):
    """Build runs the action exactly for a missing or strictly older target."""
    target_mtime = BASE + 100
    deps = [
        make_file(tmp_path / f"dep{i}", target_mtime + off)
        for i, off in enumerate(offsets)
    ]
    target = tmp_path / "target"
    if target_present:
        make_file(target, target_mtime)
    recorder = Recorder()
    task = recorder.task(target, depends=deps)
    stale = is_stale(
        target_mtime if target_present else None,
        [target_mtime + off for off in offsets],
    )

    assert task.is_up_to_date() is (not stale)
    task.build()
    assert recorder.calls == (["target"] if stale else [])
    assert task.is_up_to_date()
    task.build()
    assert recorder.calls == (["target"] if stale else [])


@pytest.mark.parametrize(
    "target_present", [False, True], ids=["target-missing", "target-present"]
)
@pytest.mark.parametrize("kinds", [("file", "task"), ("task", "task")])
@pytest.mark.parametrize(
    "offsets",
    list(itertools.product(OFFSETS, repeat=2)),
    ids=lambda combo: "deps[{:+d},{:+d}]".format(*combo),
)
def test_task_dependencies_count_like_file_dependencies(
    tmp_path, target_present, kinds, offsets
):
    """A dependency that is itself an up-to-date Task is judged by its file."""
    target_mtime = BASE + 100
    recorder = Recorder()
    deps = []
    for i, (kind, off) in enumerate(zip(kinds, offsets, strict=True)):
        path = make_file(tmp_path / f"dep{i}", target_mtime + off)
        deps.append(recorder.task(path) if kind == "task" else path)
    target = tmp_path / "target"
    if target_present:
        make_file(target, target_mtime)
    task = recorder.task(target, depends=deps)
    stale = is_stale(
        target_mtime if target_present else None,
        [target_mtime + off for off in offsets],
    )

    assert task.is_up_to_date() is (not stale)
    task.build()
    # A dependency Task that is already up to date is never rebuilt.
    assert recorder.calls == (["target"] if stale else [])
    assert task.is_up_to_date()


@pytest.mark.parametrize(
    "target_present", [False, True], ids=["target-missing", "target-present"]
)
@pytest.mark.parametrize("target_offset", [-10, 0, 10])
def test_a_missing_file_dependency_is_never_up_to_date(
    tmp_path, target_present, target_offset
):
    """Absence of any dependency makes the target stale, whatever the times."""
    present = make_file(tmp_path / "present", BASE)
    target = tmp_path / "target"
    if target_present:
        make_file(target, BASE + target_offset)
    task = Recorder().task(target, depends=[present, str(tmp_path / "absent")])

    assert not task.is_up_to_date()


@pytest.mark.parametrize("variant", ["older", "missing"])
@pytest.mark.parametrize(
    "depth,stale_level", [(d, k) for d in (2, 3, 4) for k in range(d)]
)
def test_chain_rebuilds_from_the_stale_level_upward(
    tmp_path, depth, stale_level, variant
):
    """Only the first stale level and the levels above it are rebuilt."""
    source = make_file(tmp_path / "source", BASE)
    names = [f"level{i}" for i in range(depth)]
    fresh = [BASE + 10 * (i + 1) for i in range(depth)]
    recorder = Recorder()
    tasks = []
    upstream = source
    for i, name in enumerate(names):
        path = tmp_path / name
        if i != stale_level or variant == "older":
            # Older than its dependency at the stale level, newer elsewhere.
            make_file(path, fresh[i] - 15 if i == stale_level else fresh[i])
        tasks.append(recorder.task(path, depends=upstream))
        upstream = tasks[-1]

    # Level i is up to date exactly when no level at or below i is stale.
    assert [t.is_up_to_date() for t in tasks] == [i < stale_level for i in range(depth)]
    tasks[-1].build()
    assert recorder.calls == names[stale_level:]
    assert all(t.is_up_to_date() for t in tasks)
    tasks[-1].build()
    assert recorder.calls == names[stale_level:]


@pytest.mark.parametrize("depth", [2, 3, 4])
def test_fresh_chain_runs_no_action(tmp_path, depth):
    """A chain whose every level is newer than its dependency is left alone."""
    upstream = make_file(tmp_path / "source", BASE)
    recorder = Recorder()
    task = None
    for i in range(depth):
        path = make_file(tmp_path / f"level{i}", BASE + 10 * (i + 1))
        task = recorder.task(path, depends=upstream)
        upstream = task

    assert task.is_up_to_date()
    task.build()
    assert recorder.calls == []


def test_single_string_or_task_dependency_is_not_iterated(tmp_path):
    """A bare string or Task given as depends is one dependency, not a sequence."""
    source = make_file(tmp_path / "source", BASE + 10)
    recorder = Recorder()
    from_string = recorder.task(make_file(tmp_path / "a", BASE), depends=source)
    from_task = recorder.task(make_file(tmp_path / "b", BASE + 20), depends=from_string)

    assert not from_string.is_up_to_date()
    from_task.build()
    assert recorder.calls == ["a", "b"]


@pytest.mark.parametrize(
    "target_present", [False, True], ids=["target-missing", "target-present"]
)
def test_failing_action_propagates_and_leaves_target_stale(tmp_path, target_present):
    """The action's own exception reaches the caller; the task stays stale."""
    source = make_file(tmp_path / "source", BASE + 10)
    target = tmp_path / "target"
    if target_present:
        make_file(target, BASE)
    task = Task(str(target), action=fail, depends=source)

    with pytest.raises(Boom):
        task.build()
    assert not task.is_up_to_date()


def test_failing_dependency_stops_the_build_above_it(tmp_path):
    """When a lower action fails, no higher action runs and nothing is current."""
    source = make_file(tmp_path / "source", BASE + 10)
    recorder = Recorder()
    middle = Task(str(tmp_path / "middle"), action=fail, depends=source)
    top = recorder.task(make_file(tmp_path / "top", BASE + 100), depends=middle)

    with pytest.raises(Boom):
        top.build()
    assert recorder.calls == []
    assert not middle.is_up_to_date()
    assert not top.is_up_to_date()


def test_clean_removes_only_listed_existing_files_recursively(tmp_path):
    """Listed files are removed through the whole dependency tree, nothing else."""
    source = make_file(tmp_path / "source", BASE)
    leaf_junk = make_file(tmp_path / "leaf.tmp", BASE)
    mid_junk = make_file(tmp_path / "mid.tmp", BASE)
    top_junk = make_file(tmp_path / "top.tmp", BASE)
    kept = [
        make_file(tmp_path / name, BASE)
        for name in ("unlisted.txt", "leaf", "mid", "top")
    ]
    directory = tmp_path / "listed-directory"
    directory.mkdir()
    absent = str(tmp_path / "never-created.tmp")

    recorder = Recorder()
    leaf = recorder.task(tmp_path / "leaf", depends=source, cleans=[leaf_junk, absent])
    mid = recorder.task(tmp_path / "mid", depends=[leaf, source], cleans=mid_junk)
    top = recorder.task(
        tmp_path / "top", depends=mid, cleans=[top_junk, str(directory)]
    )

    top.clean()

    assert not any(os.path.exists(p) for p in (leaf_junk, mid_junk, top_junk))
    assert all(os.path.exists(p) for p in kept + [source])
    assert directory.is_dir()
    # Cleaning again, with nothing left to remove, is not an error.
    top.clean()
    assert all(os.path.exists(p) for p in kept + [source])
    assert recorder.calls == []
