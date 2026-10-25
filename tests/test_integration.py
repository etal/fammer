"""End-to-end tests of ``fammer build``, which runs external programs.

Every test here carries the ``integration`` mark and is skipped unless the
programs it needs are on ``PATH`` and the tests run from a Git checkout, which
an sdist is not. ``pixi run test-integration`` provides the programs, runs
these tests alone and sets ``FAMMER_REQUIRE_TOOLS``, which turns each of those
skips into a failure; ``pixi run test`` leaves them out.
"""

import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
DATA = HERE / "data"
RECORDED = DATA / "py27" / "e2e"

pytestmark = pytest.mark.integration


def recorded(name):
    """One command-line run recorded from the Python 2 code base."""
    with open(RECORDED / name, encoding="utf-8") as handle:
        return json.load(handle)


def skip_unless_required(reason):
    """Skip the test, or fail it if FAMMER_REQUIRE_TOOLS is non-empty."""
    if os.environ.get("FAMMER_REQUIRE_TOOLS"):
        pytest.fail(reason)
    pytest.skip(reason)


def require(*programs):
    """Skip the test unless every program is on PATH; see skip_unless_required."""
    missing = [name for name in programs if shutil.which(name) is None]
    if missing:
        skip_unless_required("not on PATH: " + ", ".join(missing))


def run_fammer(workspace, argv):
    """Run the fammer command in workspace and return the completed process.

    ``build.align_profiles`` passes MAFFT the option ``--globalgenafpair``,
    which no MAFFT released since 2014 accepts, so the wrapper in
    ``tests/data/tool-shims`` goes first on ``PATH``, as in the reference
    capture. Remove it here once ``build`` stops passing that option.
    """
    env = dict(os.environ)
    env["PATH"] = str(DATA / "tool-shims") + os.pathsep + env["PATH"]
    return subprocess.run(["fammer", *argv], cwd=workspace, env=env,
                          capture_output=True, text=True, check=False)


def files_under(root):
    """Paths of the files under root, relative to it."""
    return {path.relative_to(root).as_posix()
            for path in root.rglob("*") if path.is_file()}


@pytest.fixture
def workspace(tmp_path):
    """A copy of the tracked files of ``profiles/SomeKinases``, as captured.

    The build writes its results beside its inputs and skips a target newer
    than its inputs, so a result left by an earlier build would change the run.
    """
    if not (REPO / ".git").exists():
        skip_unless_required("not a Git checkout, as in an sdist")
    # Only stdout is piped, so git's own error message reaches pytest's output.
    listing = subprocess.run(
        ["git", "ls-files", "-z", "--", "profiles/SomeKinases"], cwd=REPO,
        stdout=subprocess.PIPE, text=True, check=True)
    names = [name for name in listing.stdout.split("\0") if name]
    # pytest.fail, not assert: the structures test would report an
    # AssertionError raised here as its expected failure.
    if not names:
        pytest.fail("git lists no files under profiles/SomeKinases")
    for name in names:
        dest = tmp_path / Path(name).relative_to("profiles")
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(REPO / name, dest)
    return tmp_path


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
    hashed = {name: digest for name, _size, digest, kind
              in reference["manifest"] if kind == "sha256"}
    assert hashed
    differing = [name for name in sorted(hashed)
                 if hashlib.sha256((workspace / name).read_bytes()).hexdigest()
                 != hashed[name]]
    assert differing == []


# Only the assertion on the build's exit status is expected to fail; a missing
# program or checkout under FAMMER_REQUIRE_TOOLS is reported as a failure or a
# setup error.
@pytest.mark.xfail(strict=True, raises=AssertionError,
                   reason="fammer.tmalign cannot read the output of TM-align "
                   "20240303 and puts unhashable SeqRecords in its graph")
def test_build_hmmer_with_structures(workspace):
    """The same build with the structure alignments enabled."""
    require("fammer", "mafft", "hmmbuild", "hmmpress", "TMalign")
    argv = [arg for arg in recorded("build.json")["argv"] if arg != "--no-pdb"]
    proc = run_fammer(workspace, argv)
    assert proc.returncode == 0, proc.stderr
