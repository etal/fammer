# Fammer

Tools for building and curating hierarchical protein-superfamily sequence profiles.
The 2013 Python 2 code base has been ported to Python 3.12+ and current Biopython; `README.rst` is the user and developer documentation.

## Layout

- `fammer/`: the package; `cli.py` defines the `fammer` and `fammer-tmalign` entry points.
- `profiles/`: example profile trees (`SomeKinases`, `ROPK`) used as test inputs; data, not code.
- `test/`: example sequence files used as test inputs.
- `tests/`: pytest suite; `tests/data/` holds the Python 2.7 reference outputs and is documented in `tests/data/README.md`.

## Commands

Use pixi; it also installs the external programs fammer runs (MAFFT, HMMER, TM-align, MUSTANG, FastTree).

```sh
pixi run test                # unit tests, no external programs
pixi run -e py312 test       # oldest supported Python
pixi run test-integration    # build, scan, update-fasta and cluster with the external programs
pixi run lint                # ruff check
pixi run fmt                 # ruff format
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
- bioconda installs FastTree as `FastTree`, while fammer calls `fasttree`; on Linux, add a symbolic link `fasttree` in the environment's `bin`, or the `cluster` integration test fails.
- biofrills and biocma are pinned to git commits in `[tool.uv.sources]`; after changing a pin or a dependency, update both `uv.lock` and `pixi.lock`.
- The flat layout needs `packages = ["fammer"]` in `pyproject.toml`; setuptools discovery would otherwise find the example data under `profiles/` and refuse to build.

## Conventions

- Format with `ruff format`; list each formatting-only commit in `.git-blame-ignore-revs`.
- Commit messages: an imperative subject line, then a short paragraph on what changed and why; American English.
- In Markdown files, put each sentence on its own line.
