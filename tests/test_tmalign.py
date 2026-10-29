"""Properties of fammer.tmalign that need no external program.

The file-name checks of align_structs, the errors of the TM-align output
parser, and the choice of seed pairs by mst_pairs. External programs are
replaced where align_structs would start one, so that the checks are shown to
run before any alignment.
"""

import io
import itertools

import networkx
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from reference_cases import TOOL_OUTPUTS

from fammer import tmalign
from fammer._share import UsageError

# The whole standard output of TM-align 20240303 for an input file it cannot
# parse (missing, empty, or without atoms); it still exits with status 0.
UNPARSEABLE_OUTPUT = "#Total CPU time is  0.00 seconds\n"

# File names on both sides of the whitespace rule, including characters that
# a test for the space character alone would miss. The first test checks this
# classification against a FASTA round trip.
SURVIVING_NAMES = [
    "profiles/ROPK/1ATP.pdb",
    "_tmp_mustang.afasta",
    "zero\u200bwidth.pdb",
]
REJECTED_NAMES = [
    "a b.pdb",
    "dir x/s.pdb",
    "tab\there.pdb",
    "line\nbreak.pdb",
    "nbsp\xa0here.pdb",
    "ideographic\u3000space.pdb",
    "trailing.pdb ",
]

# Four structures with tied scores, so that several spanning trees share the
# maximum total and only the tie-breaking decides between them.
SCORES = {
    ("A", "B"): 0.8,
    ("A", "C"): 0.8,
    ("B", "C"): 0.8,
    ("C", "D"): 0.5,
    ("B", "D"): 0.5,
    ("A", "D"): 0.2,
}
NODES = {node for edge in SCORES for node in edge}


@pytest.fixture
def no_subprocess(monkeypatch):
    """Fail the test if align_structs starts any external program."""

    def refuse(*args, **kwargs):
        raise AssertionError(f"external program started: {args}")

    monkeypatch.setattr(tmalign.subprocess, "check_output", refuse)


def survives_fasta(name):
    """Whether a record with this identifier reads back with it from FASTA."""
    handle = io.StringIO()
    SeqIO.write(SeqRecord(Seq("ACDE"), id=name, description=""), handle, "fasta")
    handle.seek(0)
    return SeqIO.read(handle, "fasta").id == name


@pytest.mark.parametrize("name", SURVIVING_NAMES + REJECTED_NAMES)
def test_file_names_are_accepted_exactly_when_they_survive_fasta(name):
    try:
        tmalign._check_file_names(["other.pdb", name], [])
    except UsageError:
        accepted = False
    else:
        accepted = True
    assert accepted == survives_fasta(name)
    assert accepted == (name in SURVIVING_NAMES)


@pytest.mark.parametrize("name", REJECTED_NAMES)
def test_align_structs_rejects_whitespace_before_aligning(no_subprocess, name):
    """Every position is checked, the last included, before TM-align runs."""
    for names in ([name, "b.pdb", "c.pdb"], ["b.pdb", "c.pdb", name]):
        with pytest.raises(UsageError, match="whitespace") as excinfo:
            tmalign.align_structs(names)
        assert repr(name) in str(excinfo.value)


def test_align_structs_rejects_whitespace_in_seed_file_names(no_subprocess):
    with pytest.raises(UsageError, match="'my seed.fa'"):
        tmalign.align_structs(["a.pdb", "b.pdb"], ["ok.fa", "my seed.fa"])


def test_align_structs_rejects_a_repeated_structure(no_subprocess):
    with pytest.raises(UsageError, match="once: 'a.pdb'$"):
        tmalign.align_structs(["a.pdb", "b.pdb", "a.pdb"])


@pytest.mark.parametrize("names", [[], ["a.pdb"]])
def test_align_structs_of_fewer_than_two_structures_is_empty(no_subprocess, names):
    """With no pair and no seed there is nothing to align, and no program runs."""
    assert tmalign.align_structs(names) == []


def test_align_structs_rejects_missing_files_before_aligning(no_subprocess, tmp_path):
    """Only the missing names are listed, a seed file among them."""
    present, missing, seed = (
        str(tmp_path / name) for name in ("a.pdb", "b.pdb", "s.fa")
    )
    (tmp_path / "a.pdb").write_text("")
    with pytest.raises(UsageError) as excinfo:
        tmalign.align_structs([present, missing], [seed])
    assert str(excinfo.value) == f"No such file: {missing!r}, {seed!r}"


def test_output_without_a_score_raises_usage_error():
    with pytest.raises(UsageError, match="No TM-score .* a.pdb, b.pdb"):
        tmalign.read_tmalign_as_seqrec_pair(UNPARSEABLE_OUTPUT, "a.pdb", "b.pdb")


@pytest.mark.parametrize("kept_after_header", [None, 0, 1, 2])
def test_output_without_a_whole_alignment_raises_usage_error(kept_after_header):
    """The block is the header line, two sequences and the marker line between."""
    lines = (TOOL_OUTPUTS / "tmalign/1ATP__1CSN.txt").read_text().splitlines(True)
    header = next(i for i, line in enumerate(lines) if line.startswith('(":"'))
    if kept_after_header is None:
        lines = lines[:header]
    else:
        lines = lines[: header + 1 + kept_after_header]
    with pytest.raises(UsageError, match="No alignment"):
        tmalign.read_tmalign_as_seqrec_pair("".join(lines), "a.pdb", "b.pdb")


def scored_pairs(order, flipped):
    """A fresh record pair per TM-align run, as align_structs makes them."""
    pairs = []
    for left, right in order:
        if (left, right) in flipped:
            left, right = right, left
        score = SCORES[tuple(sorted((left, right)))]
        pairs.append(
            (SeqRecord(Seq("ACDE"), id=left), SeqRecord(Seq("ACDE"), id=right), score)
        )
    return pairs


def edge_set(record_pairs):
    return frozenset(frozenset((left.id, right.id)) for left, right in record_pairs)


def max_total_score():
    """The largest total score of any spanning tree, by enumeration."""
    return max(
        sum(SCORES[edge] for edge in edges)
        for edges in itertools.combinations(SCORES, len(NODES) - 1)
        if networkx.is_tree(networkx.Graph(edges))
    )


def test_mst_pairs_is_a_maximum_spanning_tree_independent_of_input_order():
    results = set()
    for order in itertools.permutations(SCORES):
        for flipped in ((), order[::2]):
            pairs = scored_pairs(order, set(flipped))
            selected = tmalign.mst_pairs(pairs)
            # The records returned are those of the input pairs, not copies.
            given = {id(rec) for left, right, _score in pairs for rec in (left, right)}
            assert all(id(rec) in given for pair in selected for rec in pair)
            results.add(edge_set(selected))
    assert len(results) == 1
    (edges,) = results
    assert len(edges) == len(NODES) - 1
    assert networkx.is_tree(networkx.Graph([tuple(edge) for edge in edges]))
    total = sum(SCORES[tuple(sorted(edge))] for edge in edges)
    assert total == pytest.approx(max_total_score())
