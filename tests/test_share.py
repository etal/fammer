"""Properties of fammer._share.parse_scanned_hits and fammer._share.write_fasta.

For each target sequence the parser keeps the profile with the smallest
best-1-domain E-value (column 8). When several profiles tie, the profile on the
first such row is kept. The tests compare the parser with an independent
oracle over every ordering of a small table.

write_fasta writes each record as one header line and one unwrapped sequence
line, with the record's identifier once in the header.
"""

import io
import itertools

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from fammer._share import parse_scanned_hits, write_fasta

HEADER = [
    "#                                  --- full sequence ---- --- best 1 dom ----\n",
    "# target name  accession  query name  accession  E-value  score  bias  E-value\n",
    "#   ----------- ---------- ----------- ---------- -------- ------ ----- -------\n",
]
FOOTER = [
    "#\n",
    "# Program:         hmmsearch\n",
    "# Query file:      profiles.hmm\n",
    "# [ok]\n",
]


def tblout_line(seq, prof, dom_eval, seq_eval="1e-20", description="-"):
    """One data line with the 18 fixed columns and a free-text description."""
    return (
        f"{seq:<24} -          {prof:<12} -          {seq_eval:>9}  100.0   0.1 "
        f"{dom_eval:>9}  99.0   0.1   1.0   1   0   0   1   1   1   1 {description}\n"
    )


def oracle(rows):
    """Expected result by direct statement of the rule, from (seq, prof, eval)."""
    expected = {}
    for seq in {row[0] for row in rows}:
        mine = [row for row in rows if row[0] == seq]
        # min() returns the first of several equal minima.
        best = min(mine, key=lambda row: float(row[2]))
        expected[seq] = (best[1], float(best[2]))
    return expected


def parse(lines):
    return parse_scanned_hits(io.StringIO("".join(lines)))


# Two sequences, one tie within each sequence and one clear minimum, written in
# different notations so that ties are numeric rather than textual.
TABLE = [
    ("seqA", "prof1", "1.0e-5"),
    ("seqA", "prof2", "1e-05"),
    ("seqA", "prof3", "2e-4"),
    ("seqB", "prof1", "3e-10"),
    ("seqB", "prof2", "1e-12"),
]


def test_every_row_order_gives_the_minimum_with_first_row_winning_ties():
    """The parser matches the oracle for all 120 orders of a five-row table."""
    for order in itertools.permutations(TABLE):
        rows = list(order)
        # Descriptions with spaces, and a full-sequence E-value that varies
        # independently of the domain E-value.
        lines = [
            tblout_line(
                seq,
                prof,
                dom,
                seq_eval=f"{i + 1}e-30",
                description=f'symbol:X species:1 "A b"  tail {i}',
            )
            for i, (seq, prof, dom) in enumerate(rows)
        ]
        assert parse(lines) == oracle(rows), rows


def test_tied_rows_keep_the_first_whatever_profile_name_or_full_evalue():
    """A tie is resolved by row position, not by profile name or full E-value."""
    lines = [
        tblout_line("s", "zzz", "1e-8", seq_eval="1e-1"),
        tblout_line("s", "aaa", "1e-8", seq_eval="1e-99"),
    ]
    assert parse(lines) == {"s": ("zzz", 1e-8)}
    assert parse(lines[::-1]) == {"s": ("aaa", 1e-8)}


def test_selection_uses_the_domain_evalue_not_the_full_sequence_evalue():
    """Column 8 decides, although column 5 orders the profiles the other way."""
    lines = [
        tblout_line("s", "wholeBest", "1e-3", seq_eval="1e-50"),
        tblout_line("s", "domainBest", "1e-9", seq_eval="1e-5"),
    ]
    assert parse(lines) == {"s": ("domainBest", 1e-9)}


def test_comment_lines_are_ignored_wherever_they_occur():
    """Header, interleaved and trailing comments contribute no hits."""
    data = [tblout_line("s1", "p1", "1e-7"), tblout_line("s2", "p2", "5e-3")]
    lines = (
        HEADER
        + [data[0], "# s3 - p9 - 1e-9 1e-9 1 0 1e-300 1 0\n"]
        + [data[1]]
        + FOOTER
    )
    assert parse(lines) == {"s1": ("p1", 1e-7), "s2": ("p2", 5e-3)}


def test_the_description_does_not_affect_the_result():
    """Descriptions with spaces, numbers and tabs are skipped as free text."""
    descriptions = [
        "-",
        "one two three",
        "1e-300 2 3 4 5 6 7 8 9 10",
        "{<Apicomplexa(e)>}Chain A, Crystal Structure",
        'symbol:CDPK1 species:3702 "Arabidopsis thaliana"  AT1G18890',
        "tab\tseparated\ttext",
    ]
    rows = [
        ("gi|254575025|pdb|3HX4|A", "Testfam", "1.6e-150"),
        ("TAIR|locus:2011201", "Testfam", "2.7e-137"),
    ]
    reference = parse([tblout_line(*row) for row in rows])
    assert reference == oracle(rows)
    for description in descriptions:
        lines = [tblout_line(*row, description=description) for row in rows]
        assert parse(lines) == reference


def test_sequences_are_independent_and_input_without_hits_is_empty():
    """Rows of one sequence never displace another; no data lines, no hits."""
    assert parse(HEADER + FOOTER) == {}
    assert parse([]) == {}
    lines = [
        tblout_line("a", "p1", "1e-2"),
        tblout_line("b", "p2", "1e-9"),
        tblout_line("a", "p3", "1e-4"),
        tblout_line("b", "p1", "1e-1"),
    ]
    assert parse(lines) == {"a": ("p3", 1e-4), "b": ("p2", 1e-9)}


@pytest.mark.parametrize(
    "text, value", [("0", 0.0), ("0.0", 0.0), ("1e-200", 1e-200), ("5.5", 5.5)]
)
def test_evalue_notations_are_returned_as_floats(text, value):
    """The E-value is a float whatever its textual form; zero is a valid value."""
    hits = parse([tblout_line("s", "p", text)])
    assert hits == {"s": ("p", value)}
    assert isinstance(hits["s"][1], float)


def test_zero_evalue_beats_every_positive_one_in_any_order():
    """An E-value of zero is a minimum, not a missing value."""
    rows = [("s", "p1", "1e-300"), ("s", "p2", "0"), ("s", "p3", "1e-100")]
    for order in itertools.permutations(rows):
        assert parse([tblout_line(*row) for row in order]) == oracle(list(order))
        assert parse([tblout_line(*row) for row in order])["s"] == ("p2", 0.0)


def written(records):
    """The text write_fasta writes for records."""
    handle = io.StringIO()
    write_fasta(records, handle)
    return handle.getvalue()


def test_a_fasta_header_is_written_back_unchanged():
    """A parsed description already starts with the identifier; it is not repeated."""
    text = ">seq1 some description\nACDEFGHIKLMNPQRSTVWY\n"
    assert written(SeqIO.parse(io.StringIO(text), "fasta")) == text


def test_a_record_without_a_description_has_the_identifier_alone():
    """No trailing space after the identifier and no blank line after the header."""
    record = SeqRecord(Seq("AC-DE"), id="s2", description="")
    assert written([record]) == ">s2\nAC-DE\n"


def test_a_clustal_record_has_the_identifier_once():
    """A Clustal parser sets the description to the identifier."""
    text = "CLUSTAL W\n\nconsensus  YLVG--RVLG\nother      YLVGK-RVLG\n"
    assert written(SeqIO.parse(io.StringIO(text), "clustal")) == (
        ">consensus\nYLVG--RVLG\n>other\nYLVGK-RVLG\n"
    )


def test_a_description_without_the_identifier_follows_it():
    """Structure alignments keep the TM-score that build.align_pdbs reads back."""
    record = SeqRecord(Seq("ACDE"), id="1ATP.pdb", description="TMalign TM-score=0.5")
    assert written([record]) == ">1ATP.pdb TMalign TM-score=0.5\nACDE\n"


def test_long_sequences_are_not_wrapped():
    """A sequence longer than the 60 columns of wrapped FASTA stays on one line."""
    record = SeqRecord(Seq("ACDEFGHIKL" * 30), id="long", description="")
    assert written([record]).splitlines() == [">long", "ACDEFGHIKL" * 30]
