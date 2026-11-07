"""Properties of fammer refine, with HMMER replaced by a table of scores.

refine holds out each sequence of a profile in turn and scores it, and the
unclassified sequences, against an HMM built without it. The lowest-scoring
sequence that scores no higher than the best unclassified one is dropped and
joins the unclassified set, and the rest are tested again, until none fails or
two sequences remain. A drop within BORDERLINE_MARGIN bits of the best
unclassified score is labeled borderline and still made. If one sequence has
the best score against two sequences that would be dropped, the profile is
kept whole, and that sequence may itself be one dropped from it earlier. A
consensus record in the unclassified set does not count. Without --dry-run,
the dropped records move from the profile's FASTA file to the unclassified
one, and no other line changes.
"""

import argparse
import logging

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from fammer import refine
from fammer._share import UsageError

# Each held-out sequence's own score, against the HMM built without it.
OWN = {"a": 50.0, "b": 200.0, "c": 40.0, "d": 300.0, "e": 250.0}
# The rival of a sequence scores this against the HMM built without it.
RIVAL = 100.0


def rival(name):
    """The unclassified sequence that scores RIVAL against the HMM without name."""
    return f"u-{name}"


@pytest.fixture
def scores(monkeypatch):
    """Replace hmmbuild and hmmsearch in refine by a table of scores.

    The HMM built without a row is named after that row. A target scores its
    own entry in OWN against the HMM that holds it out, RIVAL against the HMM
    that holds out the sequence it is the rival of, and the score set in the
    returned dict of (HMM, target) pairs; any other pair is not a hit.
    """
    extra = {}

    def leave_one_out_hmms(rows, held_out, _workdir):
        for row in rows:
            if row.id in held_out:
                yield row.id, row.id

    def get_hmm_scores(hmm, seqfname):
        with open(seqfname) as handle:
            names = [rec.id for rec in SeqIO.parse(handle, "fasta")]
        result = {}
        for name in names:
            if name == hmm and name in OWN:
                result[name] = OWN[name]
            elif name == rival(hmm):
                result[name] = RIVAL
            elif (hmm, name) in extra:
                result[name] = extra[hmm, name]
        return result

    monkeypatch.setattr(refine, "leave_one_out_hmms", leave_one_out_hmms)
    monkeypatch.setattr(refine, "get_hmm_scores", get_hmm_scores)
    return extra


def write_profile(path, names):
    """A FASTA file and the Clustal alignment that fammer build would write."""
    path.with_suffix(".fasta").write_text(
        "".join(f">{name} the {name} kinase\nACDEFG\n" for name in names)
    )
    path.with_suffix(".aln").write_text(
        "CLUSTAL X (-like) multiple sequence alignment\n\n"
        + "".join(f"{name} AC-DEFG\n" for name in names)
    )
    return str(path.with_suffix(".fasta"))


def dropped(fasta, *others):
    """The identifiers that refine_profile drops, against the rival of each
    sequence of the profile and the other unclassified sequences given."""
    names = [rival(name) for name in record_names(fasta)] + list(others)
    unclassified = [SeqRecord(Seq("ACDEFG"), id=name, description="") for name in names]
    return [rec.id for rec in refine.refine_profile(fasta, unclassified)]


def test_the_lowest_failing_sequence_goes_first_until_none_fails(
    tmp_path, scores, caplog
):
    """c (40) and a (50) score below their rivals (100); c goes first, then a,
    and each log line names the unclassified sequence it was compared with."""
    fasta = write_profile(tmp_path / "P", ["a", "b", "c", "d", "e"])
    with caplog.at_level(logging.INFO):
        assert dropped(fasta) == ["c", "a"]
    assert [m for m in caplog.messages if m.startswith("* ")] == [
        "* drop c (40.0 vs. u-c 100.0 bits)",
        "* drop a (50.0 vs. u-a 100.0 bits)",
        "* keep b (200.0 vs. u-b 100.0 bits)",
    ]


def test_a_tie_for_the_best_unclassified_score_names_the_first_identifier(
    tmp_path, scores, caplog
):
    """t-c ties u-c against the HMM without c; t-c sorts first, though listed
    after u-c."""
    scores["c", "t-c"] = RIVAL
    fasta = write_profile(tmp_path / "P", ["b", "c", "d"])
    with caplog.at_level(logging.INFO):
        assert dropped(fasta, "t-c") == ["c"]
    assert "* drop c (40.0 vs. t-c 100.0 bits)" in caplog.messages


def test_a_dropped_sequence_counts_as_unclassified(tmp_path, scores):
    """Once c is dropped, it outscores e against the HMM without e."""
    scores["e", "c"] = 260.0
    fasta = write_profile(tmp_path / "P", ["a", "b", "c", "d", "e"])
    assert dropped(fasta) == ["c", "a", "e"]


def test_a_profile_is_kept_whole_when_one_sequence_outscores_two(
    tmp_path, scores, caplog
):
    """w outscores both c and a, which then rest on one comparison; neither
    is dropped, and the warning names w and withdraws the drop of c logged
    before it."""
    scores["c", "w"] = scores["a", "w"] = 120.0
    fasta = write_profile(tmp_path / "P", ["a", "b", "c", "d", "e"])
    with caplog.at_level(logging.INFO):
        assert dropped(fasta, "w") == []
    assert "* drop c (40.0 vs. w 120.0 bits)" in caplog.messages
    assert (
        "making none of the drops above: the unclassified sequence w has the "
        "best score against both c and a; check whether it belongs in a profile"
    ) in caplog.text


def test_a_dropped_sequence_that_outscores_two_keeps_the_profile_whole(
    tmp_path, scores, caplog
):
    """c, once dropped, outscores both a and e; the warning names c as a
    sequence of the profile, not as an unclassified one."""
    scores["a", "c"] = 120.0
    scores["e", "c"] = 260.0
    fasta = write_profile(tmp_path / "P", ["a", "b", "c", "d", "e"])
    with caplog.at_level(logging.WARNING):
        assert dropped(fasta) == []
    assert (
        "c, which would be dropped from it, has the best score against both "
        "a and e; check whether the three belong together"
    ) in caplog.text


def test_two_sequences_are_always_kept(tmp_path, scores, caplog):
    """Dropping stops at two, and a profile of two is not tested at all."""
    fasta = write_profile(tmp_path / "P", ["a", "c", "b"])
    with caplog.at_level(logging.WARNING):
        assert dropped(fasta) == ["c"]
    assert "Stopped refining" in caplog.text
    fasta = write_profile(tmp_path / "Q", ["a", "c"])
    assert dropped(fasta) == []


def test_a_sequence_its_held_out_hmm_misses_fails(tmp_path, scores, monkeypatch):
    """No hit ranks below any score, but fails only beside an unclassified hit."""
    monkeypatch.delitem(OWN, "b")
    fasta = write_profile(tmp_path / "P", ["b", "d", "e"])
    assert dropped(fasta) == ["b"]
    monkeypatch.setattr(refine, "get_hmm_scores", lambda hmm, seqfname: {})
    fasta = write_profile(tmp_path / "Q", ["b", "d", "e"])
    assert dropped(fasta) == []


def test_a_score_equal_to_the_best_unclassified_one_fails(
    tmp_path, scores, monkeypatch
):
    """A sequence must score higher than its rival (100), not merely as high."""
    monkeypatch.setitem(OWN, "a", 100.0)
    fasta = write_profile(tmp_path / "P", ["a", "b", "d", "e"])
    assert dropped(fasta) == ["a"]


def test_a_drop_within_the_margin_of_the_best_unclassified_is_labeled(
    tmp_path, scores, monkeypatch, caplog
):
    """a (95) is within 10 bits of its rival (100) and is labeled borderline, but
    still dropped; c (40) and b, which its held-out HMM misses, are not."""
    monkeypatch.setitem(OWN, "a", 95.0)
    monkeypatch.delitem(OWN, "b")
    fasta = write_profile(tmp_path / "P", ["a", "b", "c", "d", "e"])
    with caplog.at_level(logging.INFO):
        assert dropped(fasta) == ["b", "c", "a"]
    assert [m for m in caplog.messages if m.startswith("* drop")] == [
        "* drop b (no hit vs. u-b 100.0 bits)",
        "* drop c (40.0 vs. u-c 100.0 bits)",
        "* drop a (95.0 vs. u-a 100.0 bits) -- borderline",
    ]


def refine_args(basedir, dry_run=False, exclude=None):
    return argparse.Namespace(
        basedir=str(basedir), dry_run=dry_run, unclassified=None, exclude=exclude
    )


# The unclassified set of the tree: the rival of each sequence of P and Q.
UNIQUE = [rival(name) for name in "abcdexyz"]


@pytest.fixture
def tree(tmp_path, scores, monkeypatch):
    """A directory with two profiles, a catch-all, and its -Unique set.

    The catch-all has no alignment, so refine would fail if it tested it.
    The -Unique file lacks a final line ending.
    """
    basedir = tmp_path / "Fam"
    basedir.mkdir()
    write_profile(basedir / "P", ["a", "b", "c", "d", "e"])
    write_profile(basedir / "Q", ["x", "y", "z"])
    (basedir / "Fam-Unclassified.fasta").write_text(">w\nACDE\n")
    (basedir / "Fam-Unique.fasta").write_text(
        "".join(f">{name} something\nACDEFG\n" for name in UNIQUE).rstrip("\n")
    )
    # x, y and z score like the strongest of P
    for name in "xyz":
        monkeypatch.setitem(OWN, name, 300.0)
    return basedir


def test_a_sequence_dropped_from_one_profile_counts_against_the_next(tree, scores):
    """c, dropped from P, outscores x against the HMM of Q without x."""
    scores["x", "c"] = 400.0
    refine.cmd_refine(refine_args(tree))
    assert record_names(tree / "Q.fasta") == ["y", "z"]
    assert record_names(tree / "Fam-Unique.fasta") == [*UNIQUE, "a", "c", "x"]


def record_names(path):
    with SeqIO.parse(path, "fasta") as records:
        return [rec.id for rec in records]


def test_refine_moves_the_dropped_records_and_nothing_else(tree):
    before = (tree / "P.fasta").read_text()
    unique = (tree / "Fam-Unique.fasta").read_text()
    refine.cmd_refine(refine_args(tree))
    by_name = {rec.split()[0]: ">" + rec for rec in before.split(">")[1:]}
    assert (tree / "P.fasta").read_text() == "".join(
        by_name[name] for name in ["b", "d", "e"]
    )
    assert (tree / "Fam-Unique.fasta").read_text() == (
        unique + "\n" + by_name["a"] + by_name["c"]
    )
    assert (tree / "Q.fasta").read_text().count(">") == 3


def test_a_dry_run_and_an_excluded_profile_change_nothing(tree):
    state = {path.name: path.read_text() for path in tree.iterdir()}
    refine.cmd_refine(refine_args(tree, dry_run=True))
    refine.cmd_refine(refine_args(tree, exclude="P"))
    assert {path.name: path.read_text() for path in tree.iterdir()} == state


def test_a_profile_kept_whole_changes_no_file_and_nothing_after_it(tree, scores):
    """w would drop c and a from P; P is kept whole, so c does not join the
    unclassified set and drop x from Q."""
    with open(tree / "Fam-Unique.fasta", "a") as handle:
        handle.write("\n>w\nACDEFG\n")
    scores["c", "w"] = scores["a", "w"] = 120.0
    scores["x", "c"] = 400.0
    state = {path.name: path.read_text() for path in tree.iterdir()}
    refine.cmd_refine(refine_args(tree))
    assert {path.name: path.read_text() for path in tree.iterdir()} == state


def test_a_profile_changed_since_the_build_is_refused(tree):
    with open(tree / "Q.fasta", "a") as handle:
        handle.write(">new\nACDEFG\n")
    with pytest.raises(UsageError, match="run fammer build first"):
        refine.cmd_refine(refine_args(tree, dry_run=True))


def test_a_sequence_in_both_sets_is_refused_before_any_change(tree):
    """P sorts first and would lose c and a, but Q shares x with -Unique."""
    with open(tree / "Fam-Unique.fasta", "a") as handle:
        handle.write("\n>x\nACDEFG\n")
    state = {path.name: path.read_text() for path in tree.iterdir()}
    with pytest.raises(UsageError, match="Sequences in both .*Q.fasta"):
        refine.cmd_refine(refine_args(tree))
    assert {path.name: path.read_text() for path in tree.iterdir()} == state


def test_a_consensus_record_in_the_unclassified_set_is_ignored(tree, scores):
    """P-consensus would outscore every member of P; only c and a move."""
    with open(tree / "Fam-Unique.fasta", "a") as handle:
        handle.write("\n>P-consensus\nACDEFG\n")
    for name in "abcde":
        scores[name, "P-consensus"] = 400.0
    refine.cmd_refine(refine_args(tree))
    assert record_names(tree / "P.fasta") == ["b", "d", "e"]
    assert record_names(tree / "Fam-Unique.fasta") == [
        *UNIQUE,
        "P-consensus",
        "a",
        "c",
    ]
