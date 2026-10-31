"""The command-line entry points, called in-process with an argument list."""

import pytest

from fammer import cli


@pytest.mark.parametrize(
    "argv, message",
    [
        (["build", "no-such-dir"], "Not a directory: no-such-dir"),
        (["scan", "no-such.hmm", "seqs.fasta"], "No HMM profile database at"),
        (["add", "no-such-dir", "seqs.fasta"], "Not a directory: no-such-dir"),
        (["refine", "no-such-dir"], "Not a directory: no-such-dir"),
        (["cluster", "no-such.aln"], "No sequence file at no-such.aln"),
        (["cluster", "seqs.txt"], "Input sequences must be a Clustal alignment"),
        (["update-fasta", "no-such.aln"], "Target must be a directory or .aln"),
    ],
)
def test_input_errors_are_reported_without_a_traceback(
    tmp_path, monkeypatch, capsys, argv, message
):
    """An input error prints one line naming the program and exits with 1."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "seqs.txt").touch()
    assert cli.main(argv) == 1
    err = capsys.readouterr().err
    assert err.startswith(f"fammer: error: {message}")
    assert "Traceback" not in err


@pytest.mark.parametrize("value", ["0", "-1e-5", "nan", "1e-400", "abc"])
def test_evalue_must_be_a_positive_number(capsys, value):
    """hmmsearch rejects a cutoff that is not positive, so the parser does too."""
    # Attached, because Python before 3.14 reads '-1e-5' as an option.
    with pytest.raises(SystemExit) as excinfo:
        cli.main(["scan", "db.hmm", "seqs.fasta", f"--evalue={value}"])
    assert excinfo.value.code == 2
    assert f"not a positive number: {value!r}" in capsys.readouterr().err


@pytest.mark.parametrize(
    "option, value, message",
    [
        ("--size", "0", "not a positive whole number"),
        ("--size", "2.5", "not a positive whole number"),
        ("--confidence", "-0.1", "not a non-negative number"),
        ("--confidence", "nan", "not a non-negative number"),
        ("--confidence", "inf", "not a non-negative number"),
    ],
)
def test_cluster_thresholds_are_checked(capsys, option, value, message):
    """A cluster size below one or a support threshold that is negative or not
    finite is rejected before any file is read."""
    with pytest.raises(SystemExit) as excinfo:
        cli.main(["cluster", "seqs.aln", f"{option}={value}"])
    assert excinfo.value.code == 2
    assert f"{message}: {value!r}" in capsys.readouterr().err


def test_tmalign_output_file_is_written_only_after_alignment(
    tmp_path, monkeypatch, capsys
):
    """A rejected input leaves no output file; one structure writes an empty one."""
    monkeypatch.chdir(tmp_path)
    out = tmp_path / "out.fasta"
    for argv, message in [
        (["a.pdb", "a.pdb"], "Each structure must be given once: 'a.pdb'"),
        (["a.pdb", "b.pdb"], "No such file: 'a.pdb', 'b.pdb'"),
    ]:
        assert cli.tmalign_main([*argv, "-o", str(out)]) == 1
        err = capsys.readouterr().err
        assert err == f"fammer-tmalign: error: {message}\n"
        assert not out.exists()

    # One structure and no seed: nothing to align, so no external program runs.
    assert cli.tmalign_main(["a.pdb", "-o", str(out)]) == 0
    assert out.read_text() == ""
    assert capsys.readouterr().err == f"tmalign: Wrote {out}\n"
    for argv in (["a.pdb"], ["a.pdb", "-o", "-"]):
        assert cli.tmalign_main(argv) == 0
        captured = capsys.readouterr()
        assert captured.out == ""
        assert "Wrote" not in captured.err
