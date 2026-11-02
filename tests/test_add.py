"""Properties of fammer.add.longest_domains and fammer.add.filter_score.

hmmsearch -A names each aligned domain ``<sequence name>/<start>-<end>``, so a
sequence with several domains appears several times. longest_domains keeps,
for each sequence, its longest domain among the domains of at least the
minimum length, and the first of two equally long ones.

filter_score keeps the hits whose bit score is at least the threshold.
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from fammer.add import filter_score, longest_domains


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
