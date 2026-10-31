#!/usr/bin/env python3
"""Replay recorded reference cases against the current fammer code.

Every library case recorded under py27/library, except those of hmmer.json,
which need HMMER, the update-fasta runs of e2e/update-fasta.json and the
cluster runs of e2e/cluster.json are run again with the current sources in
this repository, and each case whose outcome differs from the recorded one is
written to differences.json beside this script.  README.md beside it explains
each difference.

The cases and the comparison are those of tests/reference_cases.py, which
tests/test_reference.py uses as well.  The script needs a Git checkout and
the optional dependency biocma, which the repository's default pixi
environment provides.  Usage, from the repository root:

    pixi run python tests/data/replay.py
"""

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.dirname(HERE)
REPO = os.path.dirname(TESTS)
OUT = os.path.join(HERE, "differences.json")

sys.path[:0] = [REPO, TESTS]

from reference_cases import (  # noqa: E402
    RUN_DIFFERENCES,
    all_cases,
    case_outcome,
    comparable,
    matches,
    missing_requirements,
)

from fammer import update_fasta  # noqa: E402

assert os.path.dirname(os.path.abspath(update_fasta.__file__)) == os.path.join(
    REPO, "fammer"
), update_fasta.__file__


def main():
    """Replay every case and rewrite differences.json."""
    found = []
    for case in all_cases():
        missing = missing_requirements(case)
        if missing:
            sys.exit(f"{case.key} requires {', '.join(missing)}")
        current = case_outcome(case)
        if not matches(current, case.recorded):
            found.append(
                {
                    "file": case.file,
                    "section": case.section,
                    "case": case.case,
                    "recorded": comparable(case.recorded),
                    "current": current,
                }
            )
    for run_differences in RUN_DIFFERENCES.values():
        found.extend(run_differences())
    text = json.dumps(
        {"differences": found}, indent=1, separators=(",", ": "), ensure_ascii=True
    )
    with open(OUT, "w", encoding="ascii", newline="\n") as handle:
        handle.write(f"{text}\n")
    print(f"{len(found)} differing cases written to {OUT}")


if __name__ == "__main__":
    main()
