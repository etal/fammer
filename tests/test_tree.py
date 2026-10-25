"""Properties of the profile-tree functions of fammer.build.

``treeify_subdirs`` renders a directory tree as a Newick string and
``taskify_subdirs`` defines the build tasks for it. Both are checked against
the layout of the directory tree itself, read here independently: a family is
a ``*.fasta`` file whose stem is not also the name of a sibling directory, and
every directory is a group.
"""

import io
import os
import shutil
import subprocess
from pathlib import Path

import pytest
from Bio import Phylo

from fammer.build import taskify_subdirs, treeify_subdirs
from fammer.tasks import Task

REPO = Path(__file__).resolve().parent.parent

TIP, GROUP = 0, 1


def copy_tracked_profiles(name, dest):
    """Copy the Git-tracked files of profiles/<name> under dest; return its path."""
    if not (REPO / ".git").exists():
        pytest.skip("not a Git checkout, as in an sdist")
    listing = subprocess.run(
        ["git", "ls-files", "-z", "--", "profiles/" + name],
        cwd=REPO,
        stdout=subprocess.PIPE,
        text=True,
        check=True,
    )
    names = [item for item in listing.stdout.split("\0") if item]
    if not names:
        pytest.fail("git lists no files under profiles/" + name)
    for item in names:
        target = dest / Path(item).relative_to("profiles")
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(REPO / item, target)
    return str(dest / name)


def write_files(root, names):
    """Create empty files, and the directories above them, below root."""
    for name in names:
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(">s\nACDE\n")


def make_synthetic_tree(base):
    """A tree with a group FASTA beside its directory and an empty directory.

    Families: alpha, beta, a.b, x, y, z.  Groups: top, grp, sub, empty.
    ``grp.fasta``, ``sub.fasta`` and ``empty.fasta`` are group files and must
    not become tips.
    """
    top = base / "top"
    write_files(
        top,
        [
            "alpha.fasta",
            "beta.fasta",
            "a.b.fasta",
            "grp.fasta",
            "empty.fasta",
            "notes.txt",
            "other.aln",
            "grp/x.fasta",
            "grp/y.fasta",
            "grp/sub.fasta",
            "grp/readme.txt",
            "grp/sub/z.fasta",
        ],
    )
    (top / "empty").mkdir()
    return str(top)


def expected_structure(path):
    """Nested (kind, name, children) read from the file system."""
    entries = sorted(os.listdir(path))
    dirs = [e for e in entries if os.path.isdir(os.path.join(path, e))]
    children = [expected_structure(os.path.join(path, d)) for d in dirs]
    children += [
        (TIP, e[: -len(".fasta")], ())
        for e in entries
        if e.endswith(".fasta")
        and e[: -len(".fasta")] not in dirs
        and os.path.isfile(os.path.join(path, e))
    ]
    return (GROUP, os.path.basename(path), tuple(sorted(children)))


def parsed_structure(clade):
    """The same nested form from a Bio.Phylo clade.

    Newick has no empty group: ``()name`` is read as ``name`` with one unnamed,
    childless child, which is dropped here.
    """
    if not clade.clades:
        return (TIP, clade.name, ())
    children = [
        parsed_structure(c) for c in clade.clades if c.name is not None or c.clades
    ]
    return (GROUP, clade.name, tuple(sorted(children)))


def parse_newick(text):
    return Phylo.read(io.StringIO(text + ";"), "newick").root


@pytest.fixture(params=["SomeKinases", "ROPK", "synthetic"])
def tree(request, tmp_path):
    """Path of a profile tree copy: two tracked example trees and a synthetic one."""
    if request.param == "synthetic":
        return make_synthetic_tree(tmp_path)
    return copy_tracked_profiles(request.param, tmp_path)


def test_newick_tips_and_nodes_match_the_directory_tree(tree):
    """Tips are the family FASTA stems and internal nodes the directory names."""
    root = parse_newick(treeify_subdirs(tree))
    assert parsed_structure(root) == expected_structure(tree)


def test_group_fasta_beside_its_directory_is_not_a_tip(tmp_path):
    """grp.fasta, sub.fasta and empty.fasta stand for directories, not families."""
    top = make_synthetic_tree(tmp_path)
    root = parse_newick(treeify_subdirs(top))
    tips = sorted(t.name for t in root.get_terminals() if t.name is not None)
    assert tips == ["a.b", "alpha", "beta", "x", "y", "z"]
    groups = sorted(c.name for c in root.find_clades(terminal=False))
    assert groups == ["empty", "grp", "sub", "top"]


def test_trailing_separator_on_the_top_directory_does_not_change_the_tree(tree):
    """The root is named by its directory with or without a trailing slash."""
    assert treeify_subdirs(tree + "/") == treeify_subdirs(tree)


def walk(task, seen=None):
    """Yield every distinct Task reachable through Task.depends."""
    seen = {} if seen is None else seen
    if id(task) in seen:
        return
    seen[id(task)] = task
    yield task
    for dep in task.depends:
        if isinstance(dep, Task):
            yield from walk(dep, seen)


def family_fastas(top):
    """Paths of the family FASTA files, from the file system alone."""
    found = []
    for dirpath, dirnames, filenames in os.walk(top):
        for name in filenames:
            stem = name[: -len(".fasta")]
            if name.endswith(".fasta") and stem not in dirnames:
                found.append(os.path.join(dirpath, name))
    return found


def group_dirs(top):
    return [dirpath for dirpath, _dirnames, _filenames in os.walk(top)]


def test_one_aln_task_per_family_fasta_and_per_directory(tree):
    """Every family and every directory has exactly one .aln task, no other does."""
    result = taskify_subdirs(tree, hmmer=True, mapgaps=False, use_pdb=False, level=0)
    tasks = list(walk(result.hmm))
    aln_targets = sorted(t.target for t in tasks if t.target.endswith(".aln"))

    family_alns = [f[: -len(".fasta")] + ".aln" for f in family_fastas(tree)]
    group_alns = [d + ".aln" for d in group_dirs(tree)]
    assert aln_targets == sorted(family_alns + group_alns)
    # No two distinct tasks claim the same target; a group FASTA in particular
    # must not produce a second task for its directory's alignment.
    assert len(set(aln_targets)) == len(aln_targets)


def test_family_aln_tasks_depend_on_their_own_fasta(tree):
    """Each family-level .aln task has its FASTA, and only that, as its input."""
    result = taskify_subdirs(tree, hmmer=True, mapgaps=False, use_pdb=False, level=0)
    by_target = {t.target: t for t in walk(result.hmm)}
    for fasta in family_fastas(tree):
        task = by_target[fasta[: -len(".fasta")] + ".aln"]
        assert task.depends == [fasta]
