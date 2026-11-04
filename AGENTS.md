# Fammer

Tools for building and curating hierarchical protein-superfamily sequence profiles.
The 2013 Python 2 code base has been ported to Python 3.12+ and current Biopython; `README.rst` is the user and developer documentation.

## Layout

- `fammer/`: the package; `cli.py` defines the `fammer` and `fammer-tmalign` entry points.
- `profiles/`: example profile trees (`SomeKinases`, `ROPK`) used as test inputs; data, not code.
- `test/`: example sequence files used as test inputs.
- `tests/`: pytest suite; `tests/data/` holds the Python 2.7 reference outputs and is documented in `tests/data/README.md`.
- `.github/workflows/ci.yml`: CI; its Python matrix follows the classifiers in `pyproject.toml`.
- `.github/dependabot.yml`: monthly version updates of the actions in `ci.yml`; it does not update the `pixi-version` pin and is not configured for the Python dependencies.

## Design intent

"Purpose and design" in `README.rst` states what fammer is meant to do, with its sources (Talevich & Kannan 2013, MAPGAPS, Kannan et al. 2007).
Check a change to what a command classifies or writes against it before designing.

- Membership is relative: a sequence belongs to the profile that scores it highest, and an absolute score rule cannot separate related subfamilies.
  MAPGAPS classifies a best match to a group or superfamily profile at that level, while `add` assigns it among the profiles that have a FASTA file; the README states both.
- Catch-all sets (`*-Unclassified`, `*-Unclass`, `*-Unique`) are curated group-level profiles.
  They hold orphans, sparsely populated lineage-specific subfamilies that the curator collapsed into them, and consensus records of the group's subfamilies.
  They are not a queue to be emptied into subfamilies by score; moving a member into a subfamily is a curation decision.
- A command that writes to curated FASTA files should report borderline cases for the curator instead of writing them, and should offer a dry run (`-n`).
  `add` withholds and reports borderline hits (`_share.BORDERLINE_MARGIN`); `refine` marks borderline drops but still makes them.
- The trees under `profiles/` are curated data as well as test fixtures, and their classifications encode decisions the code cannot see.
  Whether a sequence or lineage belongs to a subfamily is a phylogenetic question for the maintainer, not one to settle with a new threshold.
- `profiles/SomeKinases` is a deliberately incomplete subset of the kinase superfamily, which breaks the method's assumption that the widespread subfamilies have profiles.
  Use it to test behavior, not to tune or judge classification rules; a measurement on it is weak evidence and should say so.

## Commands

Use pixi; it also installs the external programs fammer runs (MAFFT, HMMER, TM-align, MUSTANG, FastTree).

```sh
pixi run test                # unit tests, no external programs
pixi run -e py312 test       # oldest supported Python
pixi run test-integration    # build, scan, update-fasta and cluster with the external programs
pixi run lint                # ruff check
pixi run fmt                 # ruff format
pixi run fmt-check           # ruff format --check
```

`uv run pytest` also works, but it skips the integration tests unless the programs are on `PATH`.

## Reference data

- `tests/data/py27/` and `tests/data/tool-outputs/` are captured output; never edit them by hand.
  They are regenerated only by `tests/data/capture_py27.py`, as `tests/data/README.md` describes.
- A change that deliberately moves a recorded outcome is recorded by regenerating `tests/data/differences.json` with `pixi run python tests/data/replay.py`, and each new entry is explained under "Deliberate differences" in `tests/data/README.md`.
  `tests/test_reference.py` fails on an unexplained difference and on an entry that no longer differs.
- `tests/data/capture_py27.py` is frozen Python 2 code; ruff excludes it, and it must not be modernized.

## Pitfalls

- pytest turns warnings into errors.
  Biopython 1.85 and later leave a file opened by name open until the iterator is closed, so use `with SeqIO.parse(path, ...) as records:` or pass an open handle.
- `tests/data/tool-shims/mafft` exists only for `capture_py27.py`, whose Python 2 sources pass `--globalgenafpair`; the current code passes `--globalpair` (see "The MAFFT wrapper" in `tests/data/README.md`).
- Several code paths visit files in file-system order (an unsorted `glob`, `os.walk`); sort where output depends on the order, and in tests compare such outputs as sets, or as multisets where duplicates matter (see "File-system enumeration order" in `tests/data/README.md`).
- Do not name a console script `tmalign`: on case-insensitive APFS it collides with the `TMalign` binary in the same environment.
- biofrills and biocma are pinned to git commits in `[tool.uv.sources]`; after changing a pin or a dependency, update both `uv.lock` and `pixi.lock`.
- The flat layout needs `packages = ["fammer"]` in `pyproject.toml`; setuptools discovery would otherwise find the example data under `profiles/` and refuse to build.

## Conventions

- Format with `ruff format`; list each formatting-only commit in `.git-blame-ignore-revs`.
- Commit messages: an imperative subject line, then a short paragraph on what changed and why; American English.
- In Markdown files, put each sentence on its own line.
