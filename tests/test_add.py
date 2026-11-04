"""Properties of the selection, assignment, thresholds and reports of add.

hmmsearch -A names each aligned domain ``<sequence name>/<start>-<end>``, so a
sequence with several domains appears several times. longest_domains keeps,
for each sequence, its longest domain among the domains of at least the
minimum length, and the first of two equally long ones.

filter_score keeps the hits whose bit score is at least the threshold.

cmd_add gives each sequence the profile with the smallest e-value among those
with a FASTA file, a catch-all as readily as a subfamily. accepting_subfamily
selects, for report_promotions, the subfamily profile of smallest e-value at
or below a directory whose threshold a sequence reaches. inside tells whether
one directory lies at or below another, however either is spelled.

weakest_member_score scores the sequences of a profile's alignment that are
in its FASTA file, held out one at a time when there are at most ten of them.
acceptance_thresholds raises it to the best score of a member of a subfamily
in another directory.

borderline_hits picks the assigned hits that clear their threshold by less
than BORDERLINE_MARGIN bits, or that do not lead the closest other profile by
that much.

report_promotions reports, and does not move, each catch-all member other
than a consensus record that a subfamily of its group would accept, with its
score against the catch-all built without it (held_out_scores).
"""

import argparse
from collections import defaultdict

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from fammer import add
from fammer._share import UsageError
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


# Thresholds in bits; Weak has none, since its HMM misses its own sequences.
THRESHOLDS = {"Sub1": 100.0, "Sub2": 100.0, "Sub3": 300.0, "Weak": None}
# The directory of each profile: H-Unique and Other are in another group than
# the rest, and Nested in a subdirectory of the group.
DIRECTORIES = defaultdict(
    lambda: "T/G", {"H-Unique": "T/H", "Other": "T/H", "Nested": "T/G/N"}
)


def accepting(rows):
    """What accepting_subfamily selects below T/G with THRESHOLDS and DIRECTORIES."""
    return add.accepting_subfamily(
        rows, "T/G", lambda name: THRESHOLDS.get(name, 0.0), DIRECTORIES
    )


@pytest.mark.parametrize("catch_all", ["G-Unclassified", "G-Unclass", "G-Unique"])
def test_the_accepting_subfamily_of_least_evalue_is_selected(catch_all):
    """Of two accepting subfamilies, the smaller e-value wins over the higher
    score and the earlier row; the subfamily below its threshold and the
    catch-all are not selected."""
    rows = [
        ("s", catch_all, 1e-60, 200.0),
        ("s", "Sub3", 1e-55, 250.0),
        ("s", "Sub1", 1e-40, 180.0),
        ("s", "Sub2", 1e-45, 150.0),
    ]
    assert accepting(rows) == ("Sub2", 1e-45, 150.0)


def test_a_score_equal_to_the_threshold_is_accepted():
    assert accepting([("s", "Sub1", 1e-40, 100.0)]) == ("Sub1", 1e-40, 100.0)


def test_no_subfamily_is_selected_unless_one_of_the_group_accepts():
    """Below a threshold, without a threshold, a catch-all of another group,
    and a subfamily of another group are not selected."""
    rows = [
        ("s", "Sub1", 1e-40, 99.9),
        ("s", "Weak", 1e-50, 190.0),
        ("s", "H-Unique", 1e-55, 195.0),
        ("s", "Other", 1e-58, 199.0),
    ]
    assert accepting(rows) is None


def test_a_subfamily_below_the_directory_is_selected():
    assert accepting([("s", "Nested", 1e-40, 150.0)]) == ("Nested", 1e-40, 150.0)


def test_a_name_that_only_contains_a_catch_all_suffix_is_a_subfamily():
    """A name that merely contains a catch-all suffix, such as ROPK-Unclass-like,
    is a subfamily."""
    rows = [("s", "G-Unclass", 1e-70, 300.0), ("s", "ROPK-Unclass-like", 1e-60, 200.0)]
    assert accepting(rows) == ("ROPK-Unclass-like", 1e-60, 200.0)


def borderline(hits, rows):
    """{sequence: description} as borderline_hits gives it with THRESHOLDS."""
    return add.borderline_hits(hits, rows, lambda name: THRESHOLDS.get(name, 0.0))


def add_args(basedir, sequence_db):
    """The arguments of cmd_add, with the defaults of its command line."""
    return argparse.Namespace(
        basedir=str(basedir),
        sequence_db=str(sequence_db),
        evalue=1e-20,
        min_length=0.8,
        report_promotions=False,
        dry_run=False,
    )


def test_a_hit_near_its_threshold_or_a_rival_is_borderline():
    """Sub1's threshold is 100 bits and the margin 10: t clears it by 5, u
    leads Sub2 by 5, and s clears both by 10 or more."""
    rows = [
        ("s", "Sub1", 1e-50, 150.0),
        ("s", "Sub2", 1e-40, 140.0),
        ("t", "Sub1", 1e-50, 105.0),
        ("u", "Sub1", 1e-50, 150.0),
        ("u", "Sub2", 1e-49, 145.0),
    ]
    hits = {seq: (prof, ev, score) for seq, prof, ev, score in rows if prof == "Sub1"}
    found = borderline(hits, rows)
    assert sorted(found) == ["t", "u"]
    assert found["t"] == "t for Sub1 (105.0 vs. 100.0 bits)"
    assert found["u"] == "u for Sub1 (150.0 vs. 100.0 bits; Sub2 145.0)"


def test_a_catch_all_and_a_subfamily_of_its_group_are_rivals():
    """The catch-all leads Sub1 by 5 bits for s, so s is borderline; Sub1's
    rejection of t below its threshold and Weak's lack of one are not
    borderline."""
    rows = [
        ("s", "G-Unclassified", 1e-60, 200.0),
        ("s", "Sub1", 1e-50, 195.0),
        ("t", "Sub1", 1e-40, 95.0),
        ("w", "Weak", 1e-40, 95.0),
    ]
    hits = {
        seq: (prof, ev, score)
        for seq, prof, ev, score in rows
        if (seq, prof) != ("s", "Sub1")
    }
    assert borderline(hits, rows) == {
        "s": "s for G-Unclassified (200.0 vs. 0.0 bits; Sub1 195.0)"
    }


def test_a_catch_all_keeps_the_hits_that_it_matches_best(tmp_path, monkeypatch, caplog):
    """G-Unclassified matches t and u best, and Sub, of the same group, would
    accept both by score. u is appended to G-Unclassified; G-Unclassified
    rejects t for its length, and t is appended nowhere."""
    group = tmp_path / "T" / "G"
    group.mkdir(parents=True)
    for name in ("Sub", "G-Unclassified"):
        (group / f"{name}.fasta").write_text(f">{name}-m\nACDEFG\n")
        (group / f"{name}.hmm").write_text("")
    (tmp_path / "T_all.hmm").write_text("")
    (tmp_path / "target.fasta").write_text(">t\nKLMNPQ\n>u\nRSTVWY\n")
    rows = [
        ("t", "G-Unclassified", 1e-60, 200.0),
        ("t", "Sub", 1e-50, 150.0),
        ("u", "G-Unclassified", 1e-60, 200.0),
        ("u", "Sub", 1e-50, 150.0),
    ]

    def accepted_hits(_hmm, hits, min_score, index, _min_length):
        return [
            index[name] for name, score in hits if score >= min_score and name != "t"
        ]

    monkeypatch.setattr(add, "get_hmm_rows", lambda *_args: rows)
    monkeypatch.setattr(add, "acceptance_thresholds", lambda *_args: lambda _: 100.0)
    monkeypatch.setattr(add, "accepted_hits", accepted_hits)
    args = add_args(tmp_path / "T", tmp_path / "target.fasta")
    with caplog.at_level("WARNING"):
        add.cmd_add(args)
    assert caplog.records == []
    assert (group / "G-Unclassified.fasta").read_text() == (
        ">G-Unclassified-m\nACDEFG\n>u\nRSTVWY\n"
    )
    assert (group / "Sub.fasta").read_text() == ">Sub-m\nACDEFG\n"


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
    assert add.weakest_member_score(fasta, "P.hmm") == threshold


def test_a_record_added_since_the_build_does_not_count(tmp_path, hmmer):
    """A FASTA record that is not a row of the alignment cannot lower it.

    Twelve sequences are scored against the profile's own HMM, which would
    give the record "new" 10 bits.
    """
    fasta = write_profile(tmp_path, 12, extra=["new"])
    assert add.weakest_member_score(fasta, "P.hmm") == 500.0


def test_a_profile_without_an_alignment_scores_its_fasta_file(tmp_path, hmmer):
    fasta = write_profile(tmp_path, 3, extra=["new"], aligned=False)
    assert add.weakest_member_score(fasta, "P.hmm") == 10.0


# The members of each profile of a tree, by path: A and B are subfamilies of
# group G, C a subfamily of group H, and G-Unclassified the catch-all of G.
# Each profile's weakest member scores 100 bits.
NOISE_TREE = {
    "G/A": ["a0", "a1"],
    "G/B": ["b0"],
    "G/G-Unclassified": ["g0"],
    "H/C": ["c0", "C-consensus"],
}


@pytest.fixture
def noise_tree(tmp_path, monkeypatch):
    """Profiles of NOISE_TREE without alignments, scored from a table.

    weakest_member_score is replaced by a table of one score per profile, and
    get_hmm_scores by a table of (HMM, sequence) scores, set in the returned
    pair of dicts; a pair not in the table is not a hit.
    """
    profiles = {}
    for path, names in NOISE_TREE.items():
        fasta = tmp_path / f"{path}.fasta"
        fasta.parent.mkdir(exist_ok=True)
        fasta.write_text("".join(f">{name}\nACDEFG\n" for name in names))
        profiles[path.split("/")[1]] = (str(fasta), path.split("/")[1])
    weakest = {"A": 100.0, "B": 100.0, "C": 100.0, "G-Unclassified": 100.0}
    scores = {}

    def get_hmm_scores(hmm, seqfname):
        with open(seqfname) as handle:
            names = [rec.id for rec in SeqIO.parse(handle, "fasta")]
        return {name: scores[hmm, name] for name in names if (hmm, name) in scores}

    monkeypatch.setattr(add, "weakest_member_score", lambda f, hmm: weakest[hmm])
    monkeypatch.setattr(add, "get_hmm_scores", get_hmm_scores)
    return profiles, weakest, scores


def test_a_member_of_another_group_raises_the_threshold(tmp_path, noise_tree):
    """c0 of group H outscores A's weakest member; the sister subfamily B, the
    catch-all member g0 and the consensus record score higher still, but do
    not count."""
    profiles, _weakest, scores = noise_tree
    scores.update(
        {
            ("A", "c0"): 150.0,
            ("A", "b0"): 400.0,
            ("A", "g0"): 450.0,
            ("A", "C-consensus"): 500.0,
        }
    )
    threshold = add.acceptance_thresholds(profiles, str(tmp_path))
    assert threshold("A") == 150.0


def test_the_threshold_applies_to_catch_alls_and_needs_a_weakest_member(
    tmp_path, noise_tree
):
    """A catch-all is raised as a subfamily is; a profile whose HMM misses its
    own sequences has no threshold, whatever the others score; a member of
    another group that scores below the weakest member changes nothing."""
    profiles, weakest, scores = noise_tree
    weakest["B"] = None
    scores.update(
        {("G-Unclassified", "c0"): 180.0, ("B", "c0"): 300.0, ("C", "a1"): 90.0}
    )
    threshold = add.acceptance_thresholds(profiles, str(tmp_path))
    assert threshold("G-Unclassified") == 180.0
    assert threshold("B") is None
    assert threshold("C") == 100.0


# The catch-all members of the promotion tree, each with the rows of its scan:
# x goes to the subfamily Sub; y-consensus would too; z is assigned to Other,
# in another group; v matches H-Unique, the catch-all of that group, best, and
# goes to Other, though Sub would accept it too; w is accepted by no subfamily;
# and s is skipped.
PROMOTION_ROWS = [
    ("x", "G-Unclassified", 1e-60, 200.0),
    ("x", "Sub", 1e-50, 150.0),
    ("y-consensus", "Sub", 1e-60, 300.0),
    ("z", "G-Unclassified", 1e-60, 200.0),
    ("z", "Other", 1e-70, 250.0),
    ("v", "H-Unique", 1e-70, 250.0),
    ("v", "Other", 1e-50, 200.0),
    ("v", "Sub", 1e-60, 150.0),
    ("w", "G-Unclassified", 1e-60, 200.0),
    ("w", "Sub", 1e-50, 50.0),
    ("s", "Sub", 1e-60, 300.0),
]


@pytest.fixture
def promotion_tree(tmp_path, monkeypatch):
    """A group G with a subfamily and a catch-all, and a group H with both.

    The scan returns the rows of PROMOTION_ROWS for the sequences searched,
    and a profile accepts every hit at or above its threshold.
    """
    for directory in ("G", "H"):
        (tmp_path / directory).mkdir()
    (tmp_path / "G" / "Sub.fasta").write_text(">m\nACDEFG\n")
    (tmp_path / "G" / "G-Unclassified.fasta").write_text(
        "".join(f">{name} the {name} kinase\nACD\nEFG\n" for name in "xzvw")
        + ">y-consensus\nACDEFG\n>s\nACDEFG\n"
    )
    (tmp_path / "H" / "Other.fasta").write_text(">o\nACDEFG\n")
    (tmp_path / "H" / "H-Unique.fasta").write_text(">h\nACDEFG\n")
    profiles = {
        name: (str(tmp_path / directory / f"{name}.fasta"), name)
        for directory, name in [
            ("G", "Sub"),
            ("G", "G-Unclassified"),
            ("H", "Other"),
            ("H", "H-Unique"),
        ]
    }

    def get_hmm_rows(_hmmdb, seqfname, _evalue):
        with open(seqfname) as handle:
            names = {rec.id for rec in SeqIO.parse(handle, "fasta")}
        return [row for row in PROMOTION_ROWS if row[0] in names]

    def accepted_hits(_hmm, hits, min_score, index, _min_length):
        return [index[name] for name, score in hits if score >= min_score]

    monkeypatch.setattr(add, "get_hmm_rows", get_hmm_rows)
    monkeypatch.setattr(add, "accepted_hits", accepted_hits)
    return profiles


def test_promotions_are_reported_for_members_of_the_group_and_nothing_moves(
    tmp_path, promotion_tree, caplog
):
    """Only x is reported: y-consensus is a consensus record, z and v are
    assigned to a subfamily of another group, w reaches no subfamily's
    threshold, and s is skipped. The catch-all has no alignment, so x has no
    held-out score."""
    tree = [tmp_path / "G", tmp_path / "H"]
    before = {path: path.read_text() for top in tree for path in top.iterdir()}
    work = tmp_path / "work"
    work.mkdir()
    with caplog.at_level("INFO"):
        add.report_promotions(
            promotion_tree,
            "T_all.hmm",
            lambda name: 100.0,
            1e-20,
            0.8,
            skip={"s"},
            workdir=str(work),
        )
    assert [r.getMessage() for r in caplog.records] == [
        "* candidate x in G-Unclassified for Sub "
        "(150.0 vs. 100.0 bits; G-Unclassified without it not scored)"
    ]
    assert {path: path.read_text() for top in tree for path in top.iterdir()} == before


def test_a_member_is_scored_without_its_row_or_against_the_built_hmm(tmp_path, hmmer):
    """m0 is a row of the alignment, so it is scored against the HMM built
    without it (200 bits, not 500); "new" was added since the build, so it is
    scored against the profile's own HMM."""
    fasta = write_profile(tmp_path, 3, extra=["new"])
    records = [SeqRecord(Seq("ACDEFG"), id=name) for name in ("m0", "new")]
    work = tmp_path / "work"
    work.mkdir()
    scores = add.held_out_scores(fasta, "P.hmm", records, str(work))
    assert scores == {"m0": 200.0, "new": 10.0}


@pytest.mark.parametrize(("count", "aligned"), [(1, True), (3, False)])
def test_a_member_without_an_hmm_built_without_it_is_not_scored(
    tmp_path, hmmer, count, aligned
):
    """The only row of an alignment, and a member of a profile without one,
    have no HMM built without them; their own HMM was built from them."""
    fasta = write_profile(tmp_path, count, aligned=aligned)
    records = [SeqRecord(Seq("ACDEFG"), id="m0")]
    work = tmp_path / "work"
    work.mkdir()
    assert add.held_out_scores(fasta, "P.hmm", records, str(work)) == {}


@pytest.mark.parametrize(
    ("directory", "top", "expected"),
    [
        ("G", "G", True),
        ("G/N", "G", True),
        ("./G", "./G", True),
        ("./G/N", "./G", True),
        ("G/./N", "G/.", True),
        (".", ".", True),
        ("G/", "G/", True),
        ("G/N2", "G/N", False),
        ("H", "G", False),
    ],
)
def test_a_directory_is_inside_only_itself_and_its_ancestors_however_spelled(
    directory, top, expected
):
    assert add.inside(directory, top) is expected


def test_two_profiles_of_one_name_are_refused(tmp_path):
    """Hits are assigned by profile name, so two profiles cannot share one."""
    for group in ("A", "B"):
        (tmp_path / "T" / group).mkdir(parents=True)
        for suffix in ("fasta", "hmm"):
            (tmp_path / "T" / group / f"X.{suffix}").write_text("")
    (tmp_path / "T_all.hmm").write_text("")
    (tmp_path / "seqs.fasta").write_text(">s\nACDEFG\n")
    args = add_args(tmp_path / "T", tmp_path / "seqs.fasta")
    with pytest.raises(UsageError, match="Two profiles are named X"):
        add.cmd_add(args)
