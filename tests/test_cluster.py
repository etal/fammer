"""Tests of ``cluster``, on the recorded runs and on captured FastTree trees.

``py27/e2e/cluster.json`` records three runs of the Python 2 code base, with the
standard output of each MAFFT and FastTree call kept under
``tool-outputs/cluster/``.
``reference_cases.replay_cluster_run`` runs the current code on the same inputs
with that output replayed, so the tree handling, the clustering and the output
files are tested on any machine.
``tests/test_reference.py`` compares each replayed run with the recording; the
tests here check properties of the replayed runs and of the functions that
choose the clusters.
"""

import argparse
import io
import math
import os
import subprocess
import sys

import pytest
from Bio import Phylo
from reference_cases import (
    CLUSTER_RUNS,
    DATA,
    FASTTREE_NAME_END,
    cluster_input_records,
    is_checkout,
    load,
    replay_cluster_run,
)

from fammer import cluster
from fammer._share import UsageError

RECORDED = load(CLUSTER_RUNS)


def read_fasta_pairs(path):
    """The header and sequence lines of a FASTA file cluster wrote."""
    lines = path.read_text(encoding="utf-8").splitlines()
    return list(zip(lines[::2], lines[1::2], strict=True))


@pytest.fixture(params=sorted(RECORDED))
def replay(request, tmp_path):
    """The recorded case name, its scratch directory, and the replayed run."""
    if not is_checkout():
        pytest.skip("requires a Git checkout for profiles/, which an sdist lacks")
    name = request.param
    current, tools = replay_cluster_run(name, tmp_path)
    return name, tmp_path, current, tools


def test_fasttree_reads_safe_unique_names(replay):
    """FastTree is called once, on one name per input sequence, and no name
    contains a character at which FastTree would end it."""
    name, _, _, tools = replay
    assert tools.calls == len(RECORDED[name]["tool_calls"])
    (names,) = tools.tree_inputs
    assert len(names) == len(cluster_input_records(name))
    assert len(set(names)) == len(names)
    assert not any(FASTTREE_NAME_END.search(n) for n in names)


def test_every_run_succeeds(replay):
    """Each run succeeds, including CDK.aln, whose identifiers contain ':'."""
    _, _, current, _ = replay
    assert current["outcome"] == {"value": None}


def test_written_records_keep_the_input_header(replay):
    """Each record is written as its input header and one line of sequence.

    The header names the identifier once, and sequences read from an alignment
    are written without gaps, in upper case.
    """
    name, tmp_path, current, _ = replay
    by_id = {rec.id: rec for rec in cluster_input_records(name)}
    for fname in current["files_written"]:
        if fname.endswith(".xml"):
            continue
        for header, seq in read_fasta_pairs(tmp_path / fname):
            rec = by_id[header[1:].split(None, 1)[0]]
            assert header == ">" + rec.description
            if name.endswith(".aln"):
                assert seq == str(rec.seq).replace("-", "").replace(".", "").upper()
            else:
                assert seq == str(rec.seq)


def test_clusters_and_unique_partition_the_input(replay):
    """Every input identifier is written to exactly one cluster or Unique file,
    and the phyloXML tree names each input identifier once."""
    name, tmp_path, current, _ = replay
    ids = [rec.id for rec in cluster_input_records(name)]
    assert len(set(ids)) == len(ids)
    members = [i for c in current["clusters"] for i in c] + current["unique"]
    assert sorted(members) == sorted(ids)
    tree = Phylo.read(tmp_path / (name + ".xml"), "phyloxml")
    assert sorted(tip.name for tip in tree.get_terminals()) == sorted(ids)


def test_clusters_are_numbered_by_size_then_members(replay):
    """Cluster file 0 is the largest, and clusters of equal size are numbered
    in the order of their sorted member identifiers.

    CDK.aln and PKL.aln each have clusters of equal size.
    """
    _, _, current, _ = replay
    clusters = current["clusters"]
    assert all(c == sorted(c) for c in clusters)
    keys = [(-len(c), c) for c in clusters]
    assert keys == sorted(keys)


# ---------------------------------------------------------------------------
# The functions that choose the clusters, on captured FastTree trees


def captured_tree(output_file):
    """A tree FastTree wrote, from tool-outputs/cluster/."""
    text = (DATA / "tool-outputs" / "cluster" / output_file).read_text("utf-8")
    return Phylo.read(io.StringIO(text), "newick")


def test_many_tips_do_not_underflow_the_threshold():
    """On the 631 tips of PKL.aln, the cut threshold is the geometric mean of
    the tip lengths, not zero, so not every internal clade is a candidate.

    The product of the tip lengths underflows to zero.
    """
    tree = cluster.prepare_tree(captured_tree("PKL.aln.0.fasttree.out"))
    lengths = [tip.branch_length for tip in tree.get_terminals()]
    assert len(lengths) == 631
    assert math.prod(lengths) == 0.0
    assert cluster.geom_mean(lengths) == pytest.approx(0.07007, abs=1e-5)
    assert len(cluster.get_clades_to_cut(tree)) == 369
    assert len(tree.get_nonterminals()) == 418


@pytest.mark.parametrize(
    ("newick", "cut_lengths"),
    [
        # The geometric mean of the tips c, d and e is 0.2.
        ("((a:0.0,b:0.0)0.9:0.05,(c:0.1,d:0.4)0.9:0.3,(e:0.2,f:0.0)0.9:0.0);", [0.3]),
        # Every tip has length zero.
        ("((a:0.0,b:0.0)0.9:0.5,(c:0.0,d:0.0)0.9:0.0);", [0.5]),
    ],
)
def test_zero_lengths_neither_set_the_threshold_nor_are_cut(newick, cut_lengths):
    """Tips of length zero, as FastTree gives identical sequences, do not make
    the threshold zero, and a branch of length zero is never cut."""
    tree = Phylo.read(io.StringIO(newick), "newick")
    tree.root.branch_length = 0.0  # as prepare_tree leaves it
    candidates = cluster.get_clades_to_cut(tree)
    assert [cl.branch_length for cl in candidates] == cut_lengths


def test_a_tree_without_supports_collapses_nothing():
    """Without support values there is no evidence that a split is weak."""
    tree = captured_tree("DYRK.fasta.1.fasttree.out")
    for clade in tree.find_clades():
        clade.confidence = None
    internal = len(tree.get_nonterminals())
    tree = cluster.prepare_tree(tree)
    assert len(tree.get_nonterminals()) == internal


def test_clusters_of_equal_size_are_numbered_by_members(tmp_path, monkeypatch):
    """The numbering does not depend on the order in which clusters were found."""
    monkeypatch.chdir(tmp_path)
    records = "".join(f">s{i}\nACDE\n" for i in range(8))
    (tmp_path / "x.fasta").write_text(records, encoding="ascii")
    tree = Phylo.read(io.StringIO("((s4,s5,s6,s7),(s0,s1,s2,s3));"), "newick")
    first, second = tree.root.clades
    clusters = {first: {"s4", "s5", "s6", "s7"}, second: {"s0", "s1", "s2", "s3"}}
    seq_idx = cluster.read_sequences("x.fasta")
    cluster.write_clusters("x.fasta", seq_idx, tree, clusters, set())
    assert (tmp_path / "x.fasta.0").read_text("ascii").startswith(">s0\n")
    assert (tmp_path / "x.fasta.1").read_text("ascii").startswith(">s4\n")


@pytest.mark.parametrize(
    ("fname", "text"),
    [
        ("x.fasta", ">a\nACDE\n>b one\nACDF\n>a two\nACDG\n>c\nAC\n>c\nAC\n"),
        ("x.aln", "CLUSTAL W\n\na   AC-D\nb   ACED\nc   ACGD\na   ACGE\n"),
    ],
)
def test_duplicate_identifiers_are_reported_before_any_program_runs(
    tmp_path, monkeypatch, fname, text
):
    """A duplicated identifier is reported by name before MAFFT, for unaligned
    input, or FastTree is run, and nothing is written."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / fname).write_text(text, encoding="ascii")

    def no_program(*args, **kwargs):
        raise AssertionError("an external program was run")

    monkeypatch.setattr(cluster.subprocess, "check_output", no_program)
    monkeypatch.setattr(cluster, "fasttree_program", no_program)
    expected = "a, c" if fname.endswith(".fasta") else "a$"
    with pytest.raises(
        UsageError, match=f"Duplicate sequence identifiers .*: {expected}"
    ):
        cluster.cmd_cluster(argparse.Namespace(sequences=fname))
    assert sorted(os.listdir(tmp_path)) == [fname]


# ---------------------------------------------------------------------------
# Finding FastTree


def make_program(directory, name):
    path = directory / name
    path.write_text("#!/bin/sh\n", encoding="ascii")
    path.chmod(0o755)
    return path


@pytest.mark.parametrize("name", cluster.FASTTREE_NAMES)
def test_fasttree_is_found_under_either_name(tmp_path, monkeypatch, name):
    program = make_program(tmp_path, name)
    monkeypatch.setenv("PATH", str(tmp_path))
    assert os.path.samefile(cluster.fasttree_program(), program)


def test_missing_fasttree_is_reported(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(UsageError, match="fasttree, FastTree"):
        cluster.fasttree_program()


def test_cluster_does_not_need_reportlab():
    """The commands import without reportlab, which fammer does not declare."""
    code = (
        "import sys; sys.modules['reportlab'] = None; import fammer.cli, fammer.cluster"
    )
    subprocess.run([sys.executable, "-c", code], check=True)
