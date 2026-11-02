"""Properties of longest_domains, filter_score and assign_hits in fammer.add.

hmmsearch -A names each aligned domain ``<sequence name>/<start>-<end>``, so a
sequence with several domains appears several times. longest_domains keeps,
for each sequence, its longest domain among the domains of at least the
minimum length, and the first of two equally long ones.

filter_score keeps the hits whose bit score is at least the threshold.

assign_hits gives each sequence the profile with the smallest e-value, unless
that profile is a catch-all and a subfamily profile accepts the sequence.
"""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
    hits = assign_hits(rows, lambda name: THRESHOLDS.get(name, 0.0))
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
    hits = assign_hits(rows, lambda name: THRESHOLDS.get(name, 0.0))
    assert hits == {"s": ("Sub2", 1e-45, 150.0)}


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
    """ROPK-Eten1 and Unclassified-like names inside a name are subfamilies."""
    rows = [
        ("s", "ROPK-Unclass-like", 1e-60, 200.0),
        ("s", "Sub1", 1e-40, 150.0),
        ("t", "ROPK-Eten1", 1e-60, 200.0),
        ("t", "Sub1", 1e-40, 150.0),
    ]
    assert assigned(rows) == {"s": "ROPK-Unclass-like", "t": "ROPK-Eten1"}
