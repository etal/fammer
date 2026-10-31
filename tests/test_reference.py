"""Compare the current code with the outputs recorded under Python 2.7.

tests/reference_cases.py enumerates the recorded cases and defines how a
current outcome is compared with a recorded one.  A case listed in
tests/data/differences.json is expected to give the outcome listed there
instead of the recorded one; tests/data/README.md explains each entry.
"""

import json
from functools import cache

import pytest
from reference_cases import (
    DATA,
    LIBRARY,
    LIBRARY_FIXTURES,
    RUN_DIFFERENCES,
    all_cases,
    case_outcome,
    comparable,
    entry_key,
    is_checkout,
    library_cases,
    matches,
    missing_requirements,
    sections,
)


@cache
def difference_entries():
    with open(DATA / "differences.json", encoding="ascii") as handle:
        return json.load(handle)["differences"]


def library_entries():
    return [
        entry for entry in difference_entries() if entry["file"] not in RUN_DIFFERENCES
    ]


@cache
def differences():
    """The expected current outcome of each library case in differences.json."""
    return {entry_key(entry): entry["current"] for entry in library_entries()}


LIBRARY_SECTIONS = [
    pytest.param(file, section, id=f"{file.removeprefix('library/')}:{section}")
    for file in LIBRARY_FIXTURES
    for section in sections(file)
]


@pytest.mark.parametrize(("file", "section"), LIBRARY_SECTIONS)
def test_library(file, section):
    cases = [case for case in library_cases(file) if case.section == section]
    missing = sorted({name for case in cases for name in missing_requirements(case)})
    if missing:
        pytest.skip("requires " + ", ".join(missing))
    expected_changes = differences()
    found = []
    for case in cases:
        actual = case_outcome(case)
        expected = expected_changes.get(case.key, case.recorded)
        if not matches(actual, expected):
            found.append((case.case, actual))
    assert not found, f"{len(found)} of {len(cases)} cases differ, first: {found[:3]}"


def test_every_library_fixture_is_compared():
    """Each library fixture has its cases enumerated, except hmmer.json."""
    recorded = {"library/" + path.name for path in LIBRARY.glob("*.json")}
    assert recorded == {*LIBRARY_FIXTURES, "library/hmmer.json"}


def test_every_recorded_section_is_compared():
    """Each library fixture section is enumerated, and nothing else is."""
    for file in LIBRARY_FIXTURES:
        enumerated = {case.section for case in library_cases(file)}
        assert enumerated == set(sections(file)), file


@pytest.mark.parametrize("file", RUN_DIFFERENCES)
def test_runs(file):
    """The replayed runs of a fixture differ from the recorded ones as listed."""
    if not is_checkout():
        pytest.skip("requires checkout")
    expected = [e for e in difference_entries() if e["file"] == file]
    assert RUN_DIFFERENCES[file]() == expected


def test_differences_name_changed_cases():
    """Each library entry of differences.json names a recorded case that
    changed, with its recorded outcome as the fixture holds it."""
    cases = {case.key: case for case in all_cases()}
    found = []
    for entry in library_entries():
        case = cases.get(entry_key(entry))
        if (
            case is None
            or entry["recorded"] != comparable(case.recorded)
            or matches(entry["current"], case.recorded)
        ):
            found.append(entry_key(entry))
    assert not found
