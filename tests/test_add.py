"""Properties of longest_domains, filter_score, assign_hits and profile_threshold.

hmmsearch -A names each aligned domain ``<sequence name>/<start>-<end>``, so a
sequence with several domains appears several times. longest_domains keeps,
for each sequence, its longest domain among the domains of at least the
minimum length, and the first of two equally long ones.

filter_score keeps the hits whose bit score is at least the threshold.

assign_hits gives each sequence the profile with the smallest e-value, unless
that profile is a catch-all and a subfamily profile accepts the sequence. It
also returns the hits it passed on, each with the catch-all's row.

profile_threshold scores the sequences of a profile's alignment that are in
its FASTA file, held out one at a time when there are at most ten of them.
"""

import argparse

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from fammer import add
from fammer._share import UsageError
from fammer.add import assign_hits, filter_score, longest_domains


def domain(name, start, end):
    """A Stockholm-style record of an aligned domain."""
    return SeqRecord(Seq("A" * (end - start + 1)), id=f"{name}/{start}-{end}")


def selected(records, min_len_aa):
    """The (sequence name, record identifier) pairs that are kept, in order."""
    return [(acc, rec.id) for acc, rec in longest_domains(records, min_len_aa)]


def test_keeps_the_longest_domain_wherever_it_occurs():
    records = [
        domain("a", 1, 50),
        domain("b", 1, 80),
        domain("a", 61, 160),
        domain("a", 171, 230),
        domain("b", 91, 120),
    ]
    assert selected(records, 10) == [("a", "a/61-160"), ("b", "b/1-80")]


def test_keeps_the_first_of_equally_long_domains():
    records = [domain("a", 1, 100), domain("a", 201, 300)]
    assert selected(records, 10) == [("a", "a/1-100")]


def test_skips_domains_shorter_than_the_minimum():
    records = [
        domain("a", 1, 200),
        domain("a", 301, 309),
        domain("b", 1, 99),
        domain("c", 1, 100),
    ]
    assert selected(records, 100) == [("a", "a/1-200"), ("c", "c/1-100")]


def test_sequence_names_may_contain_a_slash():
    records = [domain("sp/GK070", 5, 104), domain("sp/GK070", 201, 220)]
    assert selected(records, 10) == [("sp/GK070", "sp/GK070/5-104")]


def test_filter_score_keeps_hits_at_or_above_the_threshold():
    """A hit at the threshold is kept, one a tenth of a bit below it is not."""
    hits = [("at", 100.0), ("below", 99.9), ("above", 150.0)]
    assert list(filter_score(hits, 100.0)) == [hits[0], hits[2]]


# Thresholds in bits; Weak has none, since its HMM misses its own sequences.
THRESHOLDS = {"Sub1": 100.0, "Sub2": 100.0, "Sub3": 300.0, "Weak": None}


def assigned(rows):
    """{sequence: profile} as assign_hits gives it with THRESHOLDS."""
    hits, _passed_on = assign_hits(rows, lambda name: THRESHOLDS.get(name, 0.0))
    return {seq: prof for seq, (prof, _evalue, _score) in hits.items()}


def test_a_subfamily_keeps_its_hit_whatever_other_profiles_accept():
    rows = [("s", "Sub3", 1e-50, 150.0), ("s", "Sub1", 1e-40, 400.0)]
    assert assigned(rows) == {"s": "Sub3"}


@pytest.mark.parametrize("catch_all", ["G-Unclassified", "G-Unclass", "G-Unique"])
def test_a_catch_all_hit_goes_to_the_accepting_subfamily_of_least_evalue(catch_all):
    """Of two accepting subfamilies, the smaller e-value wins over the higher
    score and the earlier row; the subfamily below its threshold is passed by."""
    rows = [
        ("s", catch_all, 1e-60, 200.0),
        ("s", "Sub3", 1e-55, 250.0),
        ("s", "Sub1", 1e-40, 180.0),
        ("s", "Sub2", 1e-45, 150.0),
    ]
    hits, passed_on = assign_hits(rows, lambda name: THRESHOLDS.get(name, 0.0))
    assert hits == {"s": ("Sub2", 1e-45, 150.0)}
    assert passed_on == {"s": (catch_all, 1e-60, 200.0)}


def test_a_score_equal_to_the_threshold_is_accepted():
    rows = [("s", "G-Unclassified", 1e-60, 200.0), ("s", "Sub1", 1e-40, 100.0)]
    assert assigned(rows) == {"s": "Sub1"}


def test_a_catch_all_keeps_a_hit_that_no_subfamily_accepts():
    """Below a threshold, without a threshold, or accepted only by another
    catch-all, the hit stays with the catch-all that matched it best."""
    rows = [
        ("s", "G-Unclassified", 1e-60, 200.0),
        ("s", "Sub1", 1e-40, 99.9),
        ("s", "Weak", 1e-50, 190.0),
        ("s", "H-Unique", 1e-55, 195.0),
    ]
    assert assigned(rows) == {"s": "G-Unclassified"}


def test_a_name_that_only_contains_a_catch_all_suffix_is_a_subfamily():
    """A name that merely contains a catch-all suffix, such as ROPK-Unclass-like,
    is a subfamily."""
    rows = [
        ("s", "ROPK-Unclass-like", 1e-60, 200.0),
        ("s", "Sub1", 1e-40, 150.0),
        ("t", "ROPK-Eten1", 1e-60, 200.0),
        ("t", "Sub1", 1e-40, 150.0),
    ]
    assert assigned(rows) == {"s": "ROPK-Unclass-like", "t": "ROPK-Eten1"}


@pytest.fixture
def hmmer(monkeypatch):
    """Replace hmmbuild and hmmsearch in add by fixed scores.

    Sequence m<i> scores 500 + i against its profile's own HMM and 200 + i
    against the HMM built without its row; "new" scores 10 against the
    profile's own HMM.
    """

    def leave_one_out_hmms(rows, held_out, _workdir):
        for row in rows:
            if row.id in held_out:
                yield row.id, "without " + row.id

    def get_hmm_scores(hmm, seqfname):
        with open(seqfname) as handle:
            names = [rec.id for rec in SeqIO.parse(handle, "fasta")]
        if hmm.startswith("without "):
            held_out = hmm[len("without ") :]
            return {held_out: 200.0 + int(held_out[1:])} if held_out in names else {}
        return {
            name: 10.0 if name == "new" else 500.0 + int(name[1:]) for name in names
        }

    monkeypatch.setattr(add, "leave_one_out_hmms", leave_one_out_hmms)
    monkeypatch.setattr(add, "get_hmm_scores", get_hmm_scores)


def write_profile(tmp_path, count, extra=(), aligned=True):
    """A profile of count aligned sequences, and extra ones in its FASTA only."""
    names = [f"m{index}" for index in range(count)]
    fasta = tmp_path / "P.fasta"
    fasta.write_text("".join(f">{name}\nACDEFG\n" for name in [*names, *extra]))
    if aligned:
        (tmp_path / "P.aln").write_text(
            "CLUSTAL X (-like) multiple sequence alignment\n\n"
            + "".join(f"{name:<4} AC-DEFG\n" for name in names)
        )
    return str(fasta)


@pytest.mark.parametrize(("count", "threshold"), [(10, 200.0), (11, 500.0)])
def test_up_to_ten_sequences_are_scored_held_out(tmp_path, hmmer, count, threshold):
    fasta = write_profile(tmp_path, count)
    assert add.profile_threshold(fasta, "P.hmm") == threshold


def test_a_record_added_since_the_build_does_not_count(tmp_path, hmmer):
    """A FASTA record that is not a row of the alignment cannot lower it.

    Twelve sequences are scored against the profile's own HMM, which would
    give the record "new" 10 bits.
    """
    fasta = write_profile(tmp_path, 12, extra=["new"])
    assert add.profile_threshold(fasta, "P.hmm") == 500.0


def test_a_profile_without_an_alignment_scores_its_fasta_file(tmp_path, hmmer):
    fasta = write_profile(tmp_path, 3, extra=["new"], aligned=False)
    assert add.profile_threshold(fasta, "P.hmm") == 10.0


def test_two_profiles_of_one_name_are_refused(tmp_path):
    """Hits are assigned by profile name, so two profiles cannot share one."""
    for group in ("A", "B"):
        (tmp_path / "T" / group).mkdir(parents=True)
        for suffix in ("fasta", "hmm"):
            (tmp_path / "T" / group / f"X.{suffix}").write_text("")
    (tmp_path / "T_all.hmm").write_text("")
    (tmp_path / "seqs.fasta").write_text(">s\nACDEFG\n")
    args = argparse.Namespace(
        basedir=str(tmp_path / "T"),
        sequence_db=str(tmp_path / "seqs.fasta"),
        evalue=1e-20,
        min_length=0.8,
    )
    with pytest.raises(UsageError, match="Two profiles are named X"):
        add.cmd_add(args)
