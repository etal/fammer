"""End-to-end tests of the fammer commands that run external programs.

Every test here carries the ``integration`` mark and is skipped unless the
programs it needs are on ``PATH`` and the tests run from a Git checkout, which
an sdist is not. ``pixi run test-integration`` provides the programs, runs
these tests alone and sets ``FAMMER_REQUIRE_TOOLS``, which turns each of those
skips into a failure; ``pixi run test`` leaves them out.

A build takes about a minute, so the tests of the commands that work on a
built profile tree share one build, made once per session; each test works on
its own copy of it.
"""

import hashlib
import logging
import os
import re
import shutil
import subprocess
from collections import Counter

import pytest
from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from reference_cases import (
    REPO,
    copy_tracked,
    is_checkout,
    load,
    log_messages,
)

from fammer import _share, add, cluster, refine, tasks

pytestmark = pytest.mark.integration


def recorded(name):
    """One command-line run recorded from the Python 2 code base."""
    return load("e2e/" + name)


def skip_unless_required(reason):
    """Skip the test, or fail it if FAMMER_REQUIRE_TOOLS is non-empty."""
    if os.environ.get("FAMMER_REQUIRE_TOOLS"):
        pytest.fail(reason)
    pytest.skip(reason)


def require(*programs):
    """Skip the test unless every program is on PATH; see skip_unless_required.

    A tuple names one program under any of several names.
    """
    missing = [
        "/".join(names)
        for names in (p if isinstance(p, tuple) else (p,) for p in programs)
        if tasks.which(names) is None
    ]
    if missing:
        skip_unless_required("not on PATH: " + ", ".join(missing))


def run_fammer(workspace, argv):
    """Run the fammer command in workspace and return the completed process."""
    return subprocess.run(
        ["fammer", *argv],
        cwd=workspace,
        capture_output=True,
        text=True,
        check=False,
    )


def files_under(root):
    """Paths of the files under root, relative to it."""
    return {
        path.relative_to(root).as_posix() for path in root.rglob("*") if path.is_file()
    }


def copy_tracked_profiles(dest):
    """Copy the tracked files of ``profiles/SomeKinases`` into dest, as captured.

    The build writes its results beside its inputs and skips a target newer
    than its inputs, so a result left by an earlier build would change the run.
    """
    if not is_checkout():
        skip_unless_required("not a Git checkout, as in an sdist")
    copy_tracked("profiles/SomeKinases", dest)


@pytest.fixture
def workspace(tmp_path):
    """A copy of the tracked files of ``profiles/SomeKinases``, as captured."""
    copy_tracked_profiles(tmp_path)
    return tmp_path


@pytest.fixture(scope="session")
def built_tree(tmp_path_factory):
    """The recorded build of ``SomeKinases``, made once and never modified.

    Tests take a copy through ``built_workspace``.
    """
    require("fammer", "mafft", "hmmbuild", "hmmpress")
    root = tmp_path_factory.mktemp("built-tree")
    copy_tracked_profiles(root)
    proc = run_fammer(root, recorded("build.json")["argv"])
    if proc.returncode != 0:
        pytest.fail(f"fammer build exited with {proc.returncode}: {proc.stderr}")
    return root


@pytest.fixture
def built_workspace(built_tree, tmp_path):
    """A private copy of the built tree, with its modification times.

    ``fammer`` decides what to rebuild from modification times, which
    ``shutil.copytree`` keeps.
    """
    workspace = tmp_path / "built"
    shutil.copytree(built_tree, workspace)
    return workspace


@pytest.fixture
def scan_workspace(built_workspace):
    """The built tree with the sequence files that the recorded scans read."""
    for name in ("twenty.fasta", "TgME49.pep"):
        shutil.copy(REPO / "test" / name, built_workspace)
    return built_workspace


def test_build_hmmer_without_structures(workspace):
    """The recorded build writes the files that the Python 2 code base wrote.

    Every created file name is compared, and the SHA-256 of each file that the
    capture hashed whole. The HMMER profiles and the pressed database are
    compared by name only: the capture hashed each profile without its
    ``DATE`` line and recorded only the size of the pressed files, which embed
    the build date. Even normalized, profiles can differ between platforms in
    the last printed digit, and ``build.all_hmm`` concatenates them in the
    order ``find`` lists them.
    """
    require("fammer", "mafft", "hmmbuild", "hmmpress")
    reference = recorded("build.json")
    inputs = files_under(workspace)
    proc = run_fammer(workspace, reference["argv"])
    assert proc.returncode == 0, proc.stderr
    assert files_under(workspace) == inputs | set(reference["created"])
    hashed = {
        name: digest
        for name, _size, digest, kind in reference["manifest"]
        if kind == "sha256"
    }
    assert hashed
    differing = [
        name
        for name in sorted(hashed)
        if hashlib.sha256((workspace / name).read_bytes()).hexdigest() != hashed[name]
    ]
    assert differing == []


def test_build_hmmer_with_structures(workspace):
    """The same build with the structure alignments enabled.

    Each structure alignment that the build writes has a row for every PDB
    file in its directory. A directory with a single PDB file has nothing to
    align, so its structure alignment is empty and the alignment of the parent
    directory has a row for that structure instead. TM-align rows are named by
    path. Where MUSTANG is installed, its alignment is a further MAFFT seed,
    whose rows are named by file name, and ``purge_seqs`` keeps only the first
    of two rows with the same gapped sequence, so a structure may be present
    under either name.
    """
    require("fammer", "mafft", "hmmbuild", "hmmpress", "TMalign")
    argv = [arg for arg in recorded("build.json")["argv"] if arg != "--no-pdb"]
    proc = run_fammer(workspace, argv)
    assert proc.returncode == 0, proc.stderr
    rows = {
        pdbseq: {os.path.basename(name) for name, _seq in read_records(pdbseq, "fasta")}
        for pdbseq in workspace.rglob("*.pdb.seq")
    }
    aligned, in_parent = set(), set()
    for pdbseq, names in sorted(rows.items()):
        topdir = pdbseq.with_name(pdbseq.name.removesuffix(".pdb.seq"))
        pdbs = {path.name for path in topdir.glob("*.pdb")}
        label = pdbseq.relative_to(workspace).as_posix()
        if names:
            assert pdbs <= names, label
            aligned.add(label)
        elif pdbs:
            parent = topdir.parent.with_name(topdir.parent.name + ".pdb.seq")
            assert pdbs <= rows[parent], label
            in_parent.add(label)
    assert {"SomeKinases.pdb.seq", "SomeKinases/CMGC.pdb.seq"} <= aligned
    assert in_parent == {"SomeKinases/STE.pdb.seq"}


def test_build_of_groups_without_families(workspace):
    """A directory holding only group subdirectories is aligned from their seeds.

    The directory has no family FASTA files of its own, so its alignment holds
    exactly the rows of the two group profiles: each group's consensus and its
    families' consensus sequences. A families file left by an earlier build
    must not contribute a row.
    """
    require("fammer", "mafft")
    ste = workspace / "SomeKinases" / "STE"
    top = workspace / "Groups"
    groups = {"GrpA": ("STE7", "STE11"), "GrpB": ("YSK_STE20", "STE-Unclassified")}
    for group, families in groups.items():
        (top / group).mkdir(parents=True)
        for family in families:
            shutil.copy(ste / f"{family}.fasta", top / group)
    (workspace / "Groups.families.fa").write_text(">Stale\n" + "A" * 50 + "\n")
    proc = run_fammer(workspace, ["build", "Groups", "--no-pdb"])
    assert proc.returncode == 0, proc.stderr
    expected = sorted(
        row
        for group in groups
        for row in ungapped_alignment(top / f"{group}.cons.seq", "fasta")
    )
    assert len(expected) == 6
    assert sorted(ungapped_alignment(workspace / "Groups.aln")) == expected


# Commands that work on the built tree.

SCAN_RUNS = recorded("scan.json")

# Eight TgME49 proteins, each with a strong hit (an E-value below 1e-60) and
# each assigned a different profile by the recorded scan of all 8,320 proteins.
TGME49_SELECTION = [
    "TGME49_203010",  # Aurora
    "TGME49_218220",  # CDK
    "TGME49_233010",  # MAPK
    "TGME49_244620",  # STE11
    "TGME49_265330",  # GSK3
    "TGME49_305860",  # CDPK
    "TGME49_306480",  # SomeKinases
    "TGME49_316150",  # ULK
]


def fasta_records(text):
    """The (header line, sequence) pairs of FASTA text, as a multiset."""
    records = []
    for line in text.splitlines():
        if line.startswith(">"):
            records.append([line, ""])
        elif line.strip():
            records[-1][1] += line.strip()
    return Counter(tuple(record) for record in records)


def read_records(path, fmt):
    """The (identifier, sequence) pairs of a sequence file."""
    with open(path, encoding="utf-8") as handle:
        return [(rec.id, str(rec.seq)) for rec in SeqIO.parse(handle, fmt)]


@pytest.mark.parametrize("name", sorted(SCAN_RUNS))
def test_scan_reproduces_recorded_run(scan_workspace, name):
    """Each recorded ``scan`` run gives its recorded output on the new build.

    The summary and the table go to standard output in a sorted order, so they
    are compared as written. The files of ``-o`` and ``-O`` list their
    sequences in the iteration order of a dictionary, which Python 2 and
    Python 3 do not share, so they are compared as sets of records. The
    recorded runs also left the hmmsearch table ``<target>.tbl`` in the
    working directory; the current code writes it to a temporary file and
    removes it, so it is the one recorded file not expected here.
    """
    require("fammer", "hmmsearch")
    run = SCAN_RUNS[name]
    before = files_under(scan_workspace)
    proc = run_fammer(scan_workspace, run["argv"])
    assert proc.returncode == run["returncode"], proc.stderr
    assert proc.stdout == run["stdout"]
    assert log_messages(proc.stderr) == log_messages(run["stderr"])
    table = run["argv"][2] + ".tbl"
    recorded_files = set(run["created"]) | set(run["modified"])
    assert table in recorded_files
    assert files_under(scan_workspace) - before == recorded_files - {table}
    for fname, recorded_text in run["written_files"].items():
        written = (scan_workspace / fname).read_text(encoding="utf-8")
        assert fasta_records(written) == fasta_records(recorded_text), fname


def test_scan_of_selected_sequences_agrees_with_scan_of_all(scan_workspace):
    """A scan of a few TgME49 proteins assigns each the profile of the full scan.

    The expected rows are those of the recorded scan of the whole file. A
    smaller target set lowers every E-value by the same factor, so a strong
    hit keeps its best profile. The sequences of ``-o`` are the selected ones.
    """
    require("fammer", "hmmsearch")
    full = dict(
        line.split("\t") for line in SCAN_RUNS["TgME49.table"]["stdout"].splitlines()
    )
    expected = {name: full[name] for name in TGME49_SELECTION}
    with open(scan_workspace / "TgME49.pep", encoding="utf-8") as handle:
        selected = [rec for rec in SeqIO.parse(handle, "fasta") if rec.id in expected]
    assert len(selected) == len(expected)
    with open(scan_workspace / "selected.pep", "w", encoding="utf-8") as handle:
        SeqIO.write(selected, handle, "fasta")

    proc = run_fammer(
        scan_workspace,
        ["scan", "SomeKinases", "selected.pep", "--table", "-o", "selected.hits.fasta"],
    )
    assert proc.returncode == 0, proc.stderr
    assert dict(line.split("\t") for line in proc.stdout.splitlines()) == expected
    assert len(proc.stdout.splitlines()) == len(expected)
    with open(scan_workspace / "selected.hits.fasta", encoding="utf-8") as handle:
        hits = Counter(
            (rec.id, rec.description, str(rec.seq))
            for rec in SeqIO.parse(handle, "fasta")
        )
    assert hits == Counter((rec.id, rec.description, str(rec.seq)) for rec in selected)


def test_scan_without_hits_writes_nothing_and_succeeds(scan_workspace):
    """A sequence unrelated to every profile gives a warning and no table."""
    require("fammer", "hmmsearch")
    (scan_workspace / "poly.fasta").write_text(">poly\n" + "A" * 60 + "\n")
    proc = run_fammer(scan_workspace, ["scan", "SomeKinases", "poly.fasta"])
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout == ""
    assert "No hits in poly.fasta" in proc.stderr


@pytest.mark.parametrize("content", [b"", b"\n"], ids=["zero-bytes", "newline"])
@pytest.mark.parametrize(
    "command, returncode, message",
    [
        ("scan", 0, "No hits in empty.fasta"),
        ("add", 1, "fammer: error: No profile hits in empty.fasta"),
    ],
)
def test_an_empty_target_takes_the_no_hits_path_of_each_command(
    scan_workspace, command, returncode, message, content
):
    """A target that is empty or holds only whitespace has no hits, and no traceback.

    hmmsearch exits with an error on a file of zero bytes, even when told to
    read FASTA, and on a file of only whitespace named ``.fasta`` unless told
    to read FASTA.
    """
    require("fammer", "hmmsearch")
    (scan_workspace / "empty.fasta").write_bytes(content)
    proc = run_fammer(scan_workspace, [command, "SomeKinases", "empty.fasta"])
    assert proc.returncode == returncode, proc.stderr
    assert message in proc.stderr
    assert "Traceback" not in proc.stderr


def alignment_pairs(root):
    """The (alignment, FASTA file) paths that share a base name under root."""
    return sorted(
        (aln, aln.with_suffix(".fasta"))
        for aln in root.rglob("*.aln")
        if aln.with_suffix(".fasta").is_file()
    )


def ungapped_alignment(path, fmt="clustal"):
    """The (identifier, ungapped upper-case sequence) pairs of an alignment."""
    return [
        (name, seq.replace("-", "").replace(".", "").upper())
        for name, seq in read_records(path, fmt)
    ]


def tree_state(root):
    """The modification time and contents of every file under root."""
    return {
        path.relative_to(root).as_posix(): (path.stat().st_mtime_ns, path.read_bytes())
        for path in root.rglob("*")
        if path.is_file()
    }


def test_update_fasta_leaves_a_fresh_build_unchanged(built_workspace):
    """After a build every FASTA file holds its alignment's ungapped records.

    So ``update-fasta`` has nothing to do: no file changes, and none is
    touched.
    """
    require("fammer")
    pairs = alignment_pairs(built_workspace)
    assert {fasta for _aln, fasta in pairs} == set(built_workspace.rglob("*.fasta"))
    for aln, fasta in pairs:
        assert set(ungapped_alignment(aln)) == set(read_records(fasta, "fasta")), fasta
    before = tree_state(built_workspace)
    proc = run_fammer(built_workspace, ["update-fasta", "SomeKinases"])
    assert proc.returncode == 0, proc.stderr
    assert tree_state(built_workspace) == before


def test_update_fasta_rewrites_stale_files_once(built_workspace):
    """``update-fasta`` rebuilds a FASTA file from its alignment, then stops.

    Every FASTA file is cut down to the first record of its alignment and made
    older than the alignment. Afterwards each file holds records of the
    alignment, ungapped; no kept sequence is contained in another; and each
    dropped sequence is contained in a kept one, as the logged purge count of
    the file requires. A sequence containing X is excluded from the two
    containment checks, because an X matches any residue in one direction
    only. A second run then changes nothing.
    """
    require("fammer")
    pairs = alignment_pairs(built_workspace)
    assert pairs
    aligned = {}
    for aln, fasta in pairs:
        aligned[fasta] = ungapped_alignment(aln)
        first_id, first_seq = aligned[fasta][0]
        fasta.write_text(f">{first_id}\n{first_seq}\n", encoding="utf-8")
        stale = aln.stat().st_mtime_ns - 10**10
        os.utime(fasta, ns=(stale, stale))

    proc = run_fammer(built_workspace, ["update-fasta", "SomeKinases"])
    assert proc.returncode == 0, proc.stderr
    purged = {
        match.group(1): int(match.group(2) or 0)
        for match in re.finditer(
            r"Updated (\S+)(?: \(purged (\d+) sequences?\))?$",
            proc.stderr,
            re.MULTILINE,
        )
    }
    for _aln, fasta in pairs:
        records = aligned[fasta]
        kept = read_records(fasta, "fasta")
        kept_seqs = [seq for _name, seq in kept]
        assert set(kept) <= set(records), fasta
        assert len(kept) == len(set(kept)), fasta
        assert len(records) - len(kept) == purged.get(
            fasta.relative_to(built_workspace).as_posix(), 0
        )
        for i, seq in enumerate(kept_seqs):
            if "X" not in seq:
                assert not any(
                    seq in other for j, other in enumerate(kept_seqs) if j != i
                ), fasta
        for _name, seq in set(records) - set(kept):
            if "X" not in seq:
                assert any(seq in other for other in kept_seqs), fasta

    written = tree_state(built_workspace)
    proc = run_fammer(built_workspace, ["update-fasta", "SomeKinases"])
    assert proc.returncode == 0, proc.stderr
    assert tree_state(built_workspace) == written


# The records that the recorded run of add appended and the first run of the
# current code does not, by file. Most are rejected by the threshold, which a
# hit must reach both by the weakest member of the profile and by the members
# of the subfamilies of other groups; TGME49_255710, which the recorded run
# appended to CDPK, is reported as borderline instead. See "Acceptance
# threshold of add", "Assignment of add hits" and "Borderline hits of add" in
# tests/data/README.md.
REMOVED_RECORDS = {
    "SomeKinases/CAMK/CDPK.fasta": {"TGME49_255710"},
    "SomeKinases/CMGC/MAPK.fasta": {"TGME49_249260", "TGME49_254630"},
    "SomeKinases/STE/STE-Unclassified.fasta": {"TGME49_218400", "TGME49_272200"},
    "SomeKinases/STE/STE11.fasta": {"TGME49_209050", "TGME49_239130", "TGME49_244620"},
}

# The records that the first run of the current code appends and the recorded
# run did not, by file: TGME49_233905, an SNF1-family kinase, matches the CAMK
# group profile best, which the recorded code did not add anywhere, and
# CAMK-Unclass best among the profiles with a FASTA file.
ADDED_RECORDS = {"SomeKinases/CAMK/CAMK-Unclass.fasta": {"TGME49_233905"}}

# The hits of the first run of add that clear their threshold by less than
# _share.BORDERLINE_MARGIN bits, or that do not lead the closest other profile
# by that much.
BORDERLINE_HITS = {
    "TGME49_204280",
    "TGME49_228750",
    "TGME49_237890",
    "TGME49_255710",
}
BORDERLINE = re.compile(
    r"add: \* borderline (\S+) for \S+ \([\d.]+ vs\. [\d.]+ bits.*\)"
)


def record_id(record):
    """The identifier of a (header line, sequence) pair."""
    return record[0][1:].split()[0]


def test_add_appends_the_recorded_hits_once(scan_workspace):
    """``add`` appends the recorded hits, with the listed exceptions, and a
    second run changes nothing.

    The first run appends to each FASTA file the records that the Python 2
    code base appended, except those in REMOVED_RECORDS, each on two lines
    where the recorded run wrapped the sequence, and the records named in
    ADDED_RECORDS. The log lines are compared as a multiset, since the tree is
    walked in file-system order; the borderline hits are reported and not
    written. The second run finds each accepted hit already in the tree, so it
    writes to no file.
    """
    require("fammer", "hmmsearch")
    run = recorded("add.json")
    before = {
        path.relative_to(scan_workspace).as_posix(): path.read_text(encoding="utf-8")
        for path in (scan_workspace / "SomeKinases").rglob("*.fasta")
    }
    proc = run_fammer(scan_workspace, run["argv"])
    assert proc.returncode == run["returncode"], proc.stderr
    appended = {}
    for name, text in before.items():
        after = (scan_workspace / name).read_text(encoding="utf-8")
        assert after.startswith(text), name
        if after != text:
            appended[name] = after[len(text) :]

    expected_files = set()
    for name in run["appended"].keys() | ADDED_RECORDS.keys():
        removed = REMOVED_RECORDS.get(name, set())
        added_ids = ADDED_RECORDS.get(name, set())
        recorded_records = fasta_records(run["appended"].get(name, ""))
        kept = Counter(
            {
                record: count
                for record, count in recorded_records.items()
                if record_id(record) not in removed
            }
        )
        if not kept and not added_ids:
            continue
        expected_files.add(name)
        appended_text = appended.get(name, "")
        records = fasta_records(appended_text)
        added = Counter(
            {
                record: count
                for record, count in records.items()
                if record_id(record) in added_ids
            }
        )
        assert {record_id(record) for record in added} == added_ids, name
        assert records - added == kept, name
        assert len(appended_text.splitlines()) == 2 * records.total(), name
    assert appended.keys() == expected_files
    messages = log_messages(proc.stderr)
    reported = [match[1] for line in messages if (match := BORDERLINE.fullmatch(line))]
    assert set(reported) == BORDERLINE_HITS
    assert [line for line in messages if not BORDERLINE.fullmatch(line)] == sorted(
        [line for line in log_messages(run["stderr"]) if ": Processing " in line]
        + [
            f"add: Adding {fasta_records(text).total()} hits to profile {name}"
            for name, text in appended.items()
        ]
    )

    written = tree_state(scan_workspace)
    proc = run_fammer(scan_workspace, run["argv"])
    assert proc.returncode == 0, proc.stderr
    assert tree_state(scan_workspace) == written


def test_add_appends_one_domain_per_sequence(scan_workspace):
    """A sequence with several domains that pass ``-L`` is appended once.

    At ``-L 0.2`` two domains of TGME49_301440, a CDPK, of 196 and 62
    residues, pass the length filter for CDPK; the longer one is kept.
    """
    require("fammer", "hmmsearch")
    proc = run_fammer(scan_workspace, ["add", "SomeKinases", "TgME49.pep", "-L", "0.2"])
    assert proc.returncode == 0, proc.stderr
    for path in (scan_workspace / "SomeKinases").rglob("*.fasta"):
        names = [name for name, _seq in read_records(path, "fasta")]
        assert _share.repeated_names(names) == [], path
    fasta = scan_workspace / "SomeKinases" / "CAMK" / "CDPK.fasta"
    added = dict(read_records(fasta, "fasta"))["TGME49_301440"]
    assert len(added) == 196


@pytest.mark.parametrize("content", [">poly\n" + "A" * 60 + "\n", ""])
def test_add_skips_a_profile_that_misses_its_own_sequences(scan_workspace, content):
    """A profile whose HMM hits none of its own sequences accepts no hit.

    Without a hit of its HMM on one of its own sequences, there is no
    threshold. The ULK sequences are replaced by one that its HMM does not
    match, or by none, while the tree's HMMs still assign a TgME49 protein to
    ULK; ``add`` warns, leaves ULK.fasta alone, and goes on to the other
    profiles.
    """
    require("fammer", "hmmsearch")
    ulk = scan_workspace / "SomeKinases" / "ULK.fasta"
    ulk.write_text(content, encoding="utf-8")
    unchanged = ulk.read_bytes()
    cdpk = scan_workspace / "SomeKinases" / "CAMK" / "CDPK.fasta"
    members = len(read_records(cdpk, "fasta"))
    proc = run_fammer(scan_workspace, ["add", "SomeKinases", "TgME49.pep"])
    assert proc.returncode == 0, proc.stderr
    assert "Skipping SomeKinases/ULK.fasta: its HMM profile hits none" in proc.stderr
    assert ulk.read_bytes() == unchanged
    assert len(read_records(cdpk, "fasta")) > members


def test_add_ends_the_last_line_of_a_profile_before_appending(scan_workspace):
    """A profile whose last line has no line ending keeps its last record.

    CDPK.fasta loses its final newline, as a hand edit can leave it. ``add``
    writes one before the recorded records, so neither the last member nor the
    first appended record is joined to the other.
    """
    require("fammer", "hmmsearch")
    name = "SomeKinases/CAMK/CDPK.fasta"
    cdpk = scan_workspace / name
    text = cdpk.read_text(encoding="utf-8").rstrip("\n")
    cdpk.write_text(text, encoding="utf-8")
    members = read_records(cdpk, "fasta")
    run = recorded("add.json")
    proc = run_fammer(scan_workspace, run["argv"])
    assert proc.returncode == 0, proc.stderr
    appended = cdpk.read_text(encoding="utf-8")[len(text) :]
    assert appended.startswith("\n>")
    removed = REMOVED_RECORDS[name]
    recorded_records = fasta_records(run["appended"][name])
    assert Counter(
        {r: n for r, n in recorded_records.items() if record_id(r) not in removed}
    ) <= fasta_records(appended)
    assert read_records(cdpk, "fasta")[: len(members)] == members


def test_add_skips_a_sequence_that_another_profile_holds(scan_workspace):
    """A hit that is already in any profile of the tree is appended nowhere.

    TGME49_233010 is a hit that the first run appends to MAPK. Placed in
    CDPK.fasta first, as a run before an earlier rebuild could have left it,
    it stays the only record of that identifier in the tree, while MAPK still
    receives its other hits.
    """
    require("fammer", "hmmsearch")
    seqid = "TGME49_233010"
    with open(scan_workspace / "TgME49.pep", encoding="utf-8") as handle:
        record = next(rec for rec in SeqIO.parse(handle, "fasta") if rec.id == seqid)
    cdpk = scan_workspace / "SomeKinases" / "CAMK" / "CDPK.fasta"
    cdpk.write_text(
        cdpk.read_text(encoding="utf-8") + f">{seqid}\n{record.seq}\n", encoding="utf-8"
    )
    mapk = scan_workspace / "SomeKinases" / "CMGC" / "MAPK.fasta"
    members = len(read_records(mapk, "fasta"))
    proc = run_fammer(scan_workspace, ["add", "SomeKinases", "TgME49.pep"])
    assert proc.returncode == 0, proc.stderr
    names = [
        name
        for path in (scan_workspace / "SomeKinases").rglob("*.fasta")
        for name, _seq in read_records(path, "fasta")
    ]
    assert names.count(seqid) == 1
    assert len(read_records(mapk, "fasta")) > members


CANDIDATE = re.compile(
    r"add: \* candidate (\S+) in (\S+) for (\S+) "
    r"\(([\d.]+) vs\. [\d.]+ bits; \S+ without it ([\d.]+)\)"
)


def test_add_reports_promotions_and_moves_nothing(scan_workspace):
    """``--report-promotions`` lists the members of STE-Unclassified that STE7
    and STE11 would accept, 15 and 5, and no consensus record, the same with
    and without ``-n``, and moves none of them.

    Held out of the catch-all's HMM, each STE7 candidate still scores higher
    against STE-Unclassified than against STE7, and each STE11 candidate lower
    than against STE11. The tree is named ``./SomeKinases``, since a directory
    spelled with a ``.`` component must still count as lying within itself.
    """
    require("fammer", "hmmsearch", "hmmbuild")
    tree = scan_workspace / "SomeKinases"
    argv = ["add", "./SomeKinases", "TgME49.pep", "--report-promotions"]
    unclassified = tree / "STE" / "STE-Unclassified.fasta"
    source = unclassified.read_bytes()
    preview = run_fammer(scan_workspace, [*argv, "-n"])
    assert preview.returncode == 0, preview.stderr
    proc = run_fammer(scan_workspace, argv)
    assert proc.returncode == 0, proc.stderr
    candidates = [
        match.groups()
        for line in log_messages(proc.stderr)
        if (match := CANDIDATE.fullmatch(line))
    ]
    assert candidates == [
        match.groups()
        for line in log_messages(preview.stderr)
        if (match := CANDIDATE.fullmatch(line))
    ]
    assert Counter(dest for _seq, _source, dest, *_scores in candidates) == {
        "STE7": 15,
        "STE11": 5,
    }
    assert {source for _seq, source, *_rest in candidates} == {"STE-Unclassified"}
    assert not any(_share.is_consensus(seq) for seq, *_rest in candidates)
    for _seq, _source, dest, score, held_out in candidates:
        assert (float(held_out) > float(score)) == (dest == "STE7")
    assert unclassified.read_bytes() == source
    names = {seq for seq, *_rest in candidates}
    for dest in ("STE7", "STE11"):
        held = {
            name for name, _seq in read_records(tree / "STE" / f"{dest}.fasta", "fasta")
        }
        assert not names & held, dest


@pytest.mark.parametrize(
    ("name", "leave_one_out"),
    [("SomeKinases/CAMK/Chk1L", True), ("SomeKinases/CAMK/CAMK2", False)],
)
def test_a_small_profile_takes_its_threshold_from_held_out_scores(
    built_workspace, monkeypatch, name, leave_one_out
):
    """A profile of up to ten sequences scores each without its own row.

    Chk1L has five sequences, each of which scores higher against the HMM
    built from all five than against the HMM built from the other four, so its
    threshold falls below the lowest score against its own HMM. CAMK2 has
    eleven, so the lowest score against its own HMM is its threshold.
    """
    require("hmmsearch", "hmmbuild")
    monkeypatch.chdir(built_workspace)
    fasta, hmm = name + ".fasta", name + ".hmm"
    assert len(read_records(fasta, "fasta")) == (5 if leave_one_out else 11)
    own = min(_share.get_hmm_scores(hmm, fasta).values())
    threshold = add.weakest_member_score(fasta, hmm)
    assert (threshold < own) if leave_one_out else (threshold == own)


def test_refine_dry_run_keeps_every_sequence_and_writes_nothing(built_workspace):
    """The recorded ``refine -n`` tests ULK and Aurora, keeps all, changes no file.

    The recorded run logged one kept seed per profile, chosen in dictionary
    order; the current one logs, for each profile, the kept sequence of lowest
    held-out score against the unclassified sequence of best score.
    """
    require("fammer", "hmmsearch", "hmmbuild")
    run = recorded("refine.json")
    written = tree_state(built_workspace)
    proc = run_fammer(built_workspace, run["argv"])
    assert proc.returncode == run["returncode"], proc.stderr
    messages = log_messages(proc.stderr)
    assert [line for line in messages if ": Refining " in line] == [
        "refine: Refining SomeKinases/Aurora.fasta ...",
        "refine: Refining SomeKinases/ULK.fasta ...",
    ]
    kept = [line.split(": * keep ")[1] for line in messages if ": * keep " in line]
    assert sorted(line.split()[0] for line in kept) == ["13031_Tt", "GK043_Gl"]
    outsiders = {line.split(" vs. ")[1].split()[0] for line in kept}
    unclassified = built_workspace / "SomeKinases" / "Testfam-Unclassified.fasta"
    assert outsiders <= set(dict(read_records(unclassified, "fasta")))
    assert not any(": * drop " in line for line in messages)
    assert tree_state(built_workspace) == written


def test_refine_moves_a_member_that_an_unclassified_sequence_matches(built_workspace):
    """A copy of ULK's weakest member in the unclassified set displaces it.

    13031_Tt has the lowest held-out score of ULK's sequences, and a copy of
    it under another name scores the same, so it is dropped first; its record
    moves to the unclassified file unchanged, and the other records of ULK stay
    as they were. Until a build realigns ULK, refine refuses to run again.
    """
    require("fammer", "hmmsearch", "hmmbuild")
    tree = built_workspace / "SomeKinases"
    ulk, unclassified = tree / "ULK.fasta", tree / "Testfam-Unclassified.fasta"
    weakest = dict(read_records(ulk, "fasta"))["13031_Tt"]
    unclassified.write_text(
        unclassified.read_text(encoding="utf-8") + f">decoy\n{weakest}\n",
        encoding="utf-8",
    )
    ulk_before = fasta_records(ulk.read_text(encoding="utf-8"))
    unclassified_before = unclassified.read_text(encoding="utf-8")
    argv = ["refine", "SomeKinases", "-u", "SomeKinases/Testfam-Unclassified.fasta"]
    proc = run_fammer(built_workspace, [*argv, "-x", "Aurora"])
    assert proc.returncode == 0, proc.stderr
    drops = [line for line in proc.stderr.splitlines() if ": * drop " in line]
    assert drops[0].split(": * drop ")[1].startswith("13031_Tt ")
    moved = fasta_records(
        unclassified.read_text(encoding="utf-8")[len(unclassified_before) :]
    )
    assert "13031_Tt" in {record_id(record) for record in moved}
    assert moved.total() == len(drops)
    assert moved <= ulk_before
    assert fasta_records(ulk.read_text(encoding="utf-8")) == ulk_before - moved
    proc = run_fammer(built_workspace, [*argv, "-n"])
    assert proc.returncode == 1
    assert "ULK.aln does not hold exactly the sequences" in proc.stderr


def test_refine_drops_the_same_on_any_number_of_threads(
    built_workspace, monkeypatch, caplog
):
    """Each thread scores its rows in files of its own, and rows keep their order.

    With a copy of 13031_Tt in the unclassified set, ULK loses 13031_Tt and is
    tested again; one thread, four, and more threads than ULK has rows drop
    the same sequences and log the same scores and identifiers.
    """
    require("hmmsearch", "hmmbuild")
    tree = built_workspace / "SomeKinases"
    ulk = str(tree / "ULK.fasta")
    with open(tree / "Testfam-Unclassified.fasta") as handle:
        unclassified = [_share.as_target(rec) for rec in SeqIO.parse(handle, "fasta")]
    weakest = dict(read_records(ulk, "fasta"))["13031_Tt"]
    unclassified.append(SeqRecord(Seq(weakest), id="decoy", description=""))
    outcomes = []
    for workers in (1, 4, 64):
        monkeypatch.setattr(refine, "WORKERS", workers)
        caplog.clear()
        with caplog.at_level(logging.INFO):
            dropped = refine.refine_profile(ulk, unclassified)
        lines = [m for m in caplog.messages if m.startswith("* ")]
        outcomes.append(([rec.id for rec in dropped], lines))
    drops, _lines = outcomes[0]
    assert drops[0] == "13031_Tt"
    assert outcomes[1] == outcomes[0]
    assert outcomes[2] == outcomes[0]


@pytest.mark.parametrize("name", ["DYRK.fasta", "CDK.fasta"])
def test_cluster_partitions_the_input_sequences(tmp_path, name):
    """``cluster`` splits a sequence set into clusters and a Unique set.

    The cluster files, numbered from 0 by decreasing size, and the Unique file
    hold every input sequence once and unchanged, and the phyloXML tree has
    one tip per input sequence. Which sequences form a cluster depends on the
    FastTree version, so membership is not checked. Two identifiers of
    CDK.fasta contain ':', at which FastTree would end a name.
    """
    require("fammer", "mafft", cluster.FASTTREE_NAMES)
    source = REPO / "profiles" / "SomeKinases" / "CMGC" / name
    shutil.copy(source, tmp_path)
    proc = run_fammer(tmp_path, ["cluster", name])
    assert proc.returncode == 0, proc.stderr

    given = read_records(tmp_path / name, "fasta")
    numbered = sorted(
        (path for path in tmp_path.glob(name + ".*") if path.suffix[1:].isdigit()),
        key=lambda path: int(path.suffix[1:]),
    )
    assert [path.suffix[1:] for path in numbered] == [
        str(i) for i in range(len(numbered))
    ]
    parts = [read_records(path, "fasta") for path in numbered]
    assert all(len(part) >= 4 for part in parts)
    assert [len(part) for part in parts] == sorted(
        (len(part) for part in parts), reverse=True
    )
    unique = tmp_path / (name + ".Unique")
    if unique.exists():
        parts.append(read_records(unique, "fasta"))
    assert parts
    assert sorted(record for part in parts for record in part) == sorted(given)

    with open(tmp_path / (name + ".xml"), encoding="utf-8") as handle:
        tree = Phylo.read(handle, "phyloxml")
    assert sorted(tip.name for tip in tree.get_terminals()) == sorted(
        seqid for seqid, _seq in given
    )


# The recorded calls ran in the workspace of the recorded build, with paths
# relative to it. The E-values are the two-digit values that hmmsearch prints,
# so an exact comparison holds across builds of one HMMER version.

HMMER_CALLS = load("library/hmmer.json")


@pytest.mark.parametrize("hmm_name", sorted(HMMER_CALLS["get_hmm_evalues"]))
def test_member_evalues_reproduce_recorded_values(
    built_workspace, monkeypatch, hmm_name
):
    """The best E-value of each sequence of a family against its own profile.

    The recorded function, ``_share.get_hmm_evalues``, was removed with the
    old ``refine``, its only caller; ``get_hmm_hits`` at ``MEMBER_EVALUE``
    returns the same E-values with the profile and the bit score.
    """
    require("hmmsearch")
    monkeypatch.chdir(built_workspace)
    fasta = hmm_name[: -len(".hmm")] + ".fasta"
    hits = _share.get_hmm_hits(hmm_name, fasta, _share.MEMBER_EVALUE)
    evalues = sorted((seq, evalue) for seq, (_prof, evalue, _score) in hits.items())
    expected = HMMER_CALLS["get_hmm_evalues"][hmm_name]["value"]
    assert evalues == [tuple(pair) for pair in expected]


def test_get_hmm_len_reproduces_recorded_lengths(built_workspace, monkeypatch):
    """Every built profile, even the combined one, has its recorded length."""
    monkeypatch.chdir(built_workspace)
    built = {
        path.relative_to(built_workspace).as_posix()
        for path in built_workspace.rglob("*.hmm")
        if "_all" not in path.name
    }
    expected = dict(HMMER_CALLS["get_hmm_len"])
    assert built == set(expected)
    assert {name: add.get_hmm_len(name) for name in built} == {
        name: call["value"] for name, call in expected.items()
    }


def test_walk_profiles_pairs_every_fasta_with_its_profile(built_workspace, monkeypatch):
    """The recorded pairs are found, in whatever order the file system lists."""
    monkeypatch.chdir(built_workspace)
    expected = HMMER_CALLS["walk_profiles"]["value"]
    assert sorted(add.walk_profiles("SomeKinases")) == sorted(
        tuple(pair) for pair in expected
    )
