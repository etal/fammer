"""Replay of the recorded ``cluster.cmd_cluster`` runs, without external programs.

``py27/e2e/cluster.json`` records three runs of the Python 2 code base, with the
standard output of each MAFFT and FastTree call kept under
``tool-outputs/cluster/``.
Here ``subprocess.check_output`` is replaced by a function that checks the
command line against the recorded one and returns the recorded output, so the
tree handling, the clustering and the output files of the current code are
compared with the recording on any machine.
"""

import argparse
import logging
import os
import shutil

import pytest
from Bio import AlignIO, SeqIO
from reference_cases import DATA, FIXTURES, REPO, is_checkout, load

from fammer import cluster

RECORDED = load("e2e/cluster.json")

# Where the capture script took each input from.
INPUTS = {
    "CDK.aln": FIXTURES / "e2e" / "build" / "SomeKinases" / "CMGC" / "CDK.aln",
    "DYRK.fasta": REPO / "profiles" / "SomeKinases" / "CMGC" / "DYRK.fasta",
    "PKL.aln": REPO / "profiles" / "ROPK" / "PKL.aln",
}

TEMPORARY_FILE = "<tmpfile>"


def read_headers(path):
    """Header lines of a FASTA file, without the leading '>'."""
    with open(path, encoding="utf-8") as handle:
        return [line[1:] for line in handle.read().splitlines() if line[0:1] == ">"]


def input_ids(name):
    """Identifiers of the sequences in an input file, as read by Biopython."""
    path = INPUTS[name]
    with open(path, encoding="utf-8") as handle:
        if name.endswith(".aln"):
            return [rec.id for rec in AlignIO.read(handle, "clustal")]
        return [rec.id for rec in SeqIO.parse(handle, "fasta")]


class RecordedTools:
    """Stand-in for ``subprocess.check_output`` that replays recorded calls."""

    def __init__(self, expected):
        self.expected = expected
        self.calls = 0
        self.temporary_inputs = []

    def __call__(self, argv, *args, **kwargs):
        assert self.calls < len(self.expected), f"unrecorded call {list(argv)}"
        recorded_argv, output_file = self.expected[self.calls]
        self.calls += 1
        assert kwargs.get("text") is True
        assert len(argv) == len(recorded_argv)
        for actual, wanted in zip(argv, recorded_argv, strict=True):
            if wanted == TEMPORARY_FILE:
                # FastTree reads this file, so it must be complete when called.
                assert os.path.isabs(actual)
                with open(actual, encoding="utf-8") as handle:
                    self.temporary_inputs.append(
                        [line[1:].rstrip("\n") for line in handle if line[0:1] == ">"]
                    )
            else:
                assert actual == wanted
        return (DATA / output_file).read_bytes().decode("utf-8")


def run_case(name, tmp_path, monkeypatch, caplog):
    """Run cmd_cluster on one recorded input in tmp_path.

    Returns the KeyError it raised, or None (any other exception propagates),
    the names of the files it wrote, the log messages at level INFO and above,
    and the stand-in for the tools.
    """
    if not is_checkout():
        pytest.skip("requires a Git checkout for profiles/, which an sdist lacks")
    record = RECORDED[name]
    shutil.copy(INPUTS[name], tmp_path)
    monkeypatch.chdir(tmp_path)
    tools = RecordedTools(record["tool_calls"])
    monkeypatch.setattr(cluster.subprocess, "check_output", tools)
    before = set(os.listdir(tmp_path))
    error = None
    with caplog.at_level(logging.INFO):
        try:
            cluster.cmd_cluster(argparse.Namespace(sequences=name))
        except KeyError as exc:
            error = exc
    written = set(os.listdir(tmp_path)) - before
    messages = [
        (rec.levelname, rec.getMessage())
        for rec in caplog.records
        if rec.levelno >= logging.INFO
    ]
    return error, written, messages, tools


@pytest.fixture(params=sorted(RECORDED))
def replay(request, tmp_path, monkeypatch, caplog):
    """The recorded case name with the outcome of replaying it."""
    name = request.param
    return name, tmp_path, run_case(name, tmp_path, monkeypatch, caplog)


@pytest.fixture
def successful_replay(replay):
    """The replay of a case whose recorded run succeeded."""
    name = replay[0]
    if "error" in RECORDED[name]["result"]:
        pytest.skip("the recorded run fails; see test_outcome_matches_the_recording")
    return replay


def test_outcome_matches_the_recording(replay):
    """A run succeeds or fails with the recorded exception type.

    The CDK.aln run fails by recorded defect: Newick reserves ':', so the tip
    names containing it come back truncated and are missing from the index.
    """
    name, _, (error, _, _, _) = replay
    result = RECORDED[name]["result"]
    if "error" in result:
        assert result["error"]["type"] == "KeyError"
        assert isinstance(error, KeyError)
    else:
        assert error is None


def test_exactly_the_recorded_tool_calls_are_made(replay):
    """MAFFT runs only for unaligned input, and FastTree reads every sequence."""
    name, _, (_, _, _, tools) = replay
    assert tools.calls == len(RECORDED[name]["tool_calls"])
    (tree_input,) = tools.temporary_inputs
    assert sorted(tree_input) == sorted(input_ids(name))


def test_files_written_match_the_recording(replay):
    """The same set of files is written, including the Unique and phyloXML ones."""
    name, _, (_, written, _, _) = replay
    assert written == set(RECORDED[name]["outputs"]["files_written"])


def test_clusters_match_the_recording(replay):
    """Clusters have the recorded members, compared without their numbers.

    Clusters of equal size are numbered in the order of a dictionary keyed by
    clade objects, which Python 2 hashed by identity, so the recording lists
    the memberships ordered by size and content.
    """
    name, tmp_path, (_, written, _, _) = replay
    outputs = RECORDED[name]["outputs"]
    found = []
    unique = None
    for fname in sorted(written):
        if fname.endswith(".xml"):
            continue
        headers = read_headers(tmp_path / fname)
        if fname.endswith(".Unique"):
            unique = headers
        else:
            found.append(headers)
    found.sort(key=lambda headers: (-len(headers), headers))
    assert found == outputs["clusters"]
    assert unique == outputs["unique"]


def test_log_matches_the_recording(replay):
    """The log holds the recorded messages, in the recorded order.

    Clusters are written by decreasing size, so the sequence of file names and
    counts does not depend on how clusters of equal size are numbered.
    """
    name, _, (_, _, messages, _) = replay
    wanted = [(level, text) for level, _, text in RECORDED[name]["result"]["log"]]
    assert messages == wanted


def test_clusters_and_unique_partition_the_input(successful_replay):
    """Every input identifier is written to exactly one cluster or Unique file."""
    name, tmp_path, (_, written, _, _) = successful_replay
    members = []
    for fname in written:
        if not fname.endswith(".xml"):
            members.extend(h.split()[0] for h in read_headers(tmp_path / fname))
    ids = input_ids(name)
    assert len(set(ids)) == len(ids)
    assert sorted(members) == sorted(ids)
    assert (tmp_path / (name + ".Unique")).is_file()


def test_cluster_numbers_follow_decreasing_size(successful_replay):
    """Cluster file 0 is the largest, and sizes never increase with the number."""
    name, tmp_path, (_, written, _, _) = successful_replay
    numbered = {
        int(fname.rsplit(".", 1)[1]): len(read_headers(tmp_path / fname))
        for fname in written
        if fname.rsplit(".", 1)[1].isdigit()
    }
    assert sorted(numbered) == list(range(len(numbered)))
    sizes = [numbered[i] for i in range(len(numbered))]
    assert sizes == sorted(sizes, reverse=True)
