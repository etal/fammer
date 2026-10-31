"""Tests of ``cluster``, on the recorded runs and on captured FastTree trees.

``py27/e2e/cluster.json`` records three runs of the Python 2 code base, with the
standard output of each MAFFT and FastTree call kept under
``tool-outputs/cluster/``.
``reference_cases.replay_cluster_run`` runs the current code on the same inputs
with that output replayed, so the tree handling, the clustering and the output
files are tested on any machine.
``tests/test_reference.py`` compares each replayed run with the recording; the
tests here check properties of the replayed runs and of the functions that
choose the clusters, the options ``--size``, ``--confidence`` and ``--tree``,
the rejection of unusable input, and the logged silhouette against
scikit-learn's.
"""

import io
import itertools
import math
import os
import shutil
import subprocess
import sys

import numpy as np
import pytest
from Bio import Phylo
from reference_cases import (
    CLUSTER_INPUTS,
    CLUSTER_RUNS,
    DATA,
    FASTTREE_NAME_END,
    cluster_args,
    cluster_input_ids,
    cluster_input_records,
    is_checkout,
    load,
    named_cluster_tree,
    outcome,
    replay_cluster_run,
    written_ids,
)
from sklearn import metrics

from fammer import cluster
from fammer._share import UsageError

RECORDED = load(CLUSTER_RUNS)


def read_fasta_pairs(path):
    """The header and sequence lines of a FASTA file cluster wrote."""
    lines = path.read_text(encoding="utf-8").splitlines()
    return list(zip(lines[::2], lines[1::2], strict=True))


@pytest.fixture
def checkout():
    """Skip unless the tests run from a Git checkout, which has profiles/."""
    if not is_checkout():
        pytest.skip("requires a Git checkout for profiles/, which an sdist lacks")


@pytest.fixture(params=sorted(RECORDED))
def replay(request, tmp_path, checkout):
    """The recorded case name, its scratch directory, and the replayed run."""
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
    ids = cluster_input_ids(name)
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


def test_the_root_support_does_not_set_the_default_threshold():
    """The root is never collapsed, so its support is not averaged with the
    splits': the mean is 0.65 rather than 0.75, and the 0.7 split is kept."""
    tree = Phylo.read(io.StringIO("((a,b)0.6,(c,d)0.7)0.95;"), "newick")
    tree = cluster.prepare_tree(tree)
    assert sorted(c.confidence for c in tree.get_nonterminals()) == [0.7, 0.95]


@pytest.mark.parametrize(
    ("output_file", "confidence", "min_size", "minimum", "count", "unclustered"),
    [
        ("CDK.aln.0.fasttree.out", None, None, 7, 8, 45),
        ("CDK.aln.0.fasttree.out", 0.9, None, 7, 8, 53),
        ("CDK.aln.0.fasttree.out", None, 4, 4, 15, 24),
        ("PKL.aln.0.fasttree.out", 0.9, None, 13, 13, 289),
        ("PKL.aln.0.fasttree.out", None, 8, 8, 26, 208),
    ],
)
def test_find_clusters_partitions_the_tips(
    output_file, confidence, min_size, minimum, count, unclustered
):
    """Every tip is in exactly one cluster or left out, and no cluster is
    smaller than the minimum size, which is explicit or, by default, 7 for the
    145 tips of CDK.aln and 13 for the 631 tips of PKL.aln."""
    tree = cluster.prepare_tree(captured_tree(output_file), confidence)
    clusters, left_out = cluster.find_clusters(tree, min_size)
    members = [name for taxa in clusters.values() for name in taxa]
    assert (len(clusters), len(left_out)) == (count, unclustered)
    assert len(members) == len(set(members))
    assert not set(members) & left_out
    assert set(members) | left_out == {tip.name for tip in tree.get_terminals()}
    assert min(len(taxa) for taxa in clusters.values()) >= minimum


def test_the_default_minimum_size():
    """The default is the tip count to the power 0.4, rounded: 7 for the 145
    tips of CDK.aln, and 4 rather than 2 for a tree of 4 tips."""
    tree = cluster.prepare_tree(captured_tree("CDK.aln.0.fasttree.out"))
    assert cluster.find_clusters(tree) == cluster.find_clusters(tree, 7)
    assert cluster.find_clusters(tree) != cluster.find_clusters(tree, 8)
    newick = "((a:1,b:1)0.9:1,(c:1,d:1)0.9:1);"
    small = cluster.prepare_tree(Phylo.read(io.StringIO(newick), "newick"))
    assert cluster.find_clusters(small) == ({}, {"a", "b", "c", "d"})


@pytest.mark.parametrize(
    ("newick", "min_confidence", "kept"),
    [
        # On bootstrap percentages, 0 collapses only the split without a
        # support value, and 2 also the split of support 0.
        ("((a,b)0,((c,d)95,(e,f))88)100;", 0, [0, 88, 95]),
        ("((a,b)0,((c,d)95,(e,f))88)100;", 2, [88, 95]),
        # On FastTree's scale, 1 keeps only the splits of support 1.
        ("((a,b)0.5,((c,d)1.0,(e,f))0.8);", 1, [1.0]),
    ],
)
def test_thresholds_on_the_tree_scale_are_accepted(newick, min_confidence, kept):
    tree = Phylo.read(io.StringIO(newick), "newick")
    tree = cluster.prepare_tree(tree, min_confidence)
    splits = [c for c in tree.get_nonterminals() if c is not tree.root]
    assert sorted(c.confidence for c in splits) == kept


def test_clusters_of_equal_size_are_numbered_by_members(tmp_path, monkeypatch):
    """The numbering does not depend on the order in which clusters were found."""
    monkeypatch.chdir(tmp_path)
    records = "".join(f">s{i}\nACDE\n" for i in range(8))
    (tmp_path / "x.fasta").write_text(records, encoding="ascii")
    tree = Phylo.read(io.StringIO("((s4,s5,s6,s7),(s0,s1,s2,s3));"), "newick")
    first, second = tree.root.clades
    clusters = {first: {"s4", "s5", "s6", "s7"}, second: {"s0", "s1", "s2", "s3"}}
    seq_idx = cluster.read_sequences("x.fasta")
    ordered = cluster.ordered_clusters(clusters)
    cluster.write_clusters("x.fasta", seq_idx, tree, ordered, set())
    assert (tmp_path / "x.fasta.0").read_text("ascii").startswith(">s0\n")
    assert (tmp_path / "x.fasta.1").read_text("ascii").startswith(">s4\n")


def no_program(*args, **kwargs):
    raise AssertionError("an external program was run")


@pytest.fixture
def no_programs(monkeypatch):
    """Fail on any call of MAFFT or FastTree."""
    monkeypatch.setattr(cluster.subprocess, "check_output", no_program)
    monkeypatch.setattr(cluster, "fasttree_program", no_program)


@pytest.mark.parametrize(
    ("fname", "text"),
    [
        ("x.fasta", ">a\nACDE\n>b one\nACDF\n>a two\nACDG\n>c\nAC\n>c\nAC\n"),
        ("x.aln", "CLUSTAL W\n\na   AC-D\nb   ACED\nc   ACGD\na   ACGE\n"),
    ],
)
def test_duplicate_identifiers_are_reported_before_any_program_runs(
    tmp_path, monkeypatch, no_programs, fname, text
):
    """A duplicated identifier is reported by name before MAFFT, for unaligned
    input, or FastTree is run, and nothing is written."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / fname).write_text(text, encoding="ascii")
    expected = "a, c" if fname.endswith(".fasta") else "a$"
    with pytest.raises(
        UsageError, match=f"Duplicate sequence identifiers .*: {expected}"
    ):
        cluster.cmd_cluster(cluster_args(fname))
    assert sorted(os.listdir(tmp_path)) == [fname]


# ---------------------------------------------------------------------------
# The options and the silhouette


@pytest.fixture
def tree_run(tmp_path, monkeypatch, checkout, no_programs):
    """A function running cluster with --tree on the input of a recorded run
    and the tree FastTree built for it, in a scratch directory; it returns the
    outcome, with the log, and the directory."""
    workdir = tmp_path / "tree-run"
    workdir.mkdir()
    monkeypatch.chdir(workdir)

    def run(name, **options):
        shutil.copy(CLUSTER_INPUTS[name], workdir)
        (workdir / "tree.nwk").write_text(named_cluster_tree(name), encoding="utf-8")
        args = cluster_args(name, tree="tree.nwk", **options)
        return outcome(cluster.cmd_cluster, args, log=True), workdir

    return run


def test_tree_option_reproduces_the_built_tree_run(replay, tree_run):
    """Given the tree FastTree built, --tree writes the same files and the same
    log without running MAFFT or FastTree, for unaligned input too, and reads
    quoted identifiers containing ':' (CDK.aln)."""
    name, replay_dir, current, _ = replay
    result, workdir = tree_run(name)
    assert result["log"] == current["log"]
    for fname in current["files_written"]:
        if not fname.endswith(".xml"):
            assert (workdir / fname).read_bytes() == (replay_dir / fname).read_bytes()


@pytest.mark.parametrize("name", ["CDK.aln", "DYRK.fasta"])
def test_logged_silhouette_matches_sklearn(tmp_path, checkout, name):
    """The logged silhouette is scikit-learn's, from the patristic distances of
    the tree as FastTree built it, over the clustered sequences alone.

    Collapsing the weakly supported splits first would give 0.393 and 0.596
    instead of 0.398 and 0.630; including Unique would change the counts.
    """
    current, _ = replay_cluster_run(name, tmp_path)
    tree = Phylo.read(io.StringIO(named_cluster_tree(name)), "newick")
    tips = {tip.name: tip for tip in tree.get_terminals()}
    clusters = current["clusters"]
    members = [seqid for members in clusters for seqid in members]
    labels = np.repeat(np.arange(len(clusters)), [len(c) for c in clusters])
    distances = np.zeros((len(members), len(members)))
    for i, j in itertools.combinations(range(len(members)), 2):
        distances[i, j] = distances[j, i] = tree.distance(
            tips[members[i]], tips[members[j]]
        )
    score = metrics.silhouette_score(distances, labels, metric="precomputed")
    samples = metrics.silhouette_samples(distances, labels, metric="precomputed")
    expected = [
        f"Mean silhouette {score:.3f} over the {len(members)} sequences of "
        f"{len(clusters)} clusters, leaving out {len(current['unique'])} Unique "
        "sequences"
    ] + [
        f"Silhouette of {name}.{index}: {samples[labels == index].mean():.3f}"
        for index in range(len(clusters))
    ]
    logged = [message for _, _, message in current["log"] if "ilhouette" in message]
    assert logged == expected


@pytest.mark.parametrize(
    ("points", "labels"),
    [
        # A group of one point.
        ([0.0, 1.0, 5.0, 6.0, 20.0], [0, 0, 1, 1, 2]),
        # Every distance is zero.
        ([0.0, 0.0, 0.0, 0.0], [0, 0, 1, 1]),
        # Labels that are not 0 to k - 1, in no particular order.
        ([0.0, 1.0, 2.0, 10.0, 11.0], ["b", "b", "a", "a", "a"]),
    ],
)
def test_silhouette_conventions_match_sklearn(points, labels):
    x = np.array(points)
    distances = np.abs(x[:, None] - x)
    expected = metrics.silhouette_samples(distances, labels, metric="precomputed")
    assert cluster.silhouette(distances, labels) == pytest.approx(expected, abs=1e-12)


def test_patristic_distances_match_biopython():
    """Every pair, on a tree with a polytomy, a root branch, and a branch
    without a length."""
    newick = "((a:1,b:2)0.9:0.5,(c:1,d,e:3)0.8:2,(f:4)0.7:0.25):0.7;"
    tree = Phylo.read(io.StringIO(newick), "newick")
    names, distances = cluster.patristic_distances(tree)
    tips = tree.get_terminals()
    assert names == [tip.name for tip in tips]
    expected = [
        [tree.distance(a, b) if a is not b else 0.0 for b in tips] for a in tips
    ]
    assert distances == pytest.approx(np.array(expected), abs=1e-12)


def test_size_and_confidence_options_reach_the_command(tree_run):
    """Both options change which files are written, and the logged silhouette
    covers the sequences of those files."""
    name = "CDK.aln"
    count = 13
    result, workdir = tree_run(name, size=4, confidence=0.9)
    assert result["value"] is None
    written = os.listdir(workdir)
    assert sum(f"{name}.{index}" in written for index in range(count + 1)) == count
    sizes = [len(written_ids(workdir / f"{name}.{i}")) for i in range(count)]
    assert min(sizes) >= 4
    assert len(written_ids(workdir / f"{name}.Unique")) == 41
    assert f"over the {sum(sizes)} sequences of {count} clusters" in str(result["log"])


def test_fewer_than_two_clusters_have_no_silhouette(tree_run):
    result, workdir = tree_run("DYRK.fasta", size=20)
    assert (workdir / "DYRK.fasta.0").exists()
    assert not (workdir / "DYRK.fasta.1").exists()
    messages = [message for _, _, message in result["log"]]
    assert messages[-1] == "No silhouette: it needs at least two clusters (found 1)"
    assert not any(message.startswith("Silhouette") for message in messages)


FOUR_SEQUENCES = ">a\nACDE\n>b\nACDF\n>c:1\nACDG\n>d\nACDH\n"


@pytest.mark.parametrize(
    ("newick", "options", "message"),
    [
        ("((a,b),('c:1',e));", {}, "Tips of tree.nwk that name no input sequence: e$"),
        ("((a,b),'c:1');", {}, "Input sequences that no tip of tree.nwk names: d$"),
        ("((a,b),('c:1',d),a);", {}, "Duplicate tip names in tree.nwk: a$"),
        (
            "((a,b),(c:1:0.5,d));",
            {},
            "Tips of tree.nwk that name no input sequence: c$",
        ),
        ("((a,b),('c:1',d);", {}, "Could not read one Newick tree from tree.nwk"),
        ("", {}, "Could not read one Newick tree from tree.nwk"),
        (
            "((a:1,b:1):1,('c:1':1,d:1):1);",
            {"confidence": 0.5},
            "No split of the tree has a support value",
        ),
        (
            "((a,b)0.8,('c:1',d)0.6);",
            {"confidence": 0.9},
            "No split of the tree has support of 0.9 or more; the largest is 0.8,",
        ),
        # The root is not a split, and collapse_all leaves it in place.
        (
            "((a,b)0.6,('c:1',d)0.7)0.95;",
            {"confidence": 0.9},
            "No split of the tree has support of 0.9 or more; the largest is 0.7,",
        ),
        ("((a,b),('c:1',d),'');", {}, "Tips of tree.nwk without a name: 1$"),
        # A threshold up to 1 on a tree of bootstrap percentages.
        (
            "((a,b)95,('c:1',d)72)100;",
            {"confidence": 0.9},
            "A threshold of 0.9 would collapse almost no split of this tree, "
            "whose support values run up to 95;",
        ),
        (
            "((a,b)95,('c:1',d)72)100;",
            {"confidence": 1},
            "A threshold of 1 would collapse almost no split",
        ),
    ],
)
def test_unusable_trees_are_reported(
    tmp_path, monkeypatch, no_programs, newick, options, message
):
    """A tree that does not name exactly the input sequences, cannot be read,
    or does not support the threshold is reported, and nothing is written."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "x.fasta").write_text(FOUR_SEQUENCES, encoding="ascii")
    (tmp_path / "tree.nwk").write_text(newick, encoding="ascii")
    with pytest.raises(UsageError, match=message):
        cluster.cmd_cluster(cluster_args("x.fasta", tree="tree.nwk", **options))
    assert sorted(os.listdir(tmp_path)) == ["tree.nwk", "x.fasta"]


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
