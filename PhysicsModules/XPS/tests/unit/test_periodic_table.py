"""Characterization tests for xps_neo.periodic_table.

periodic_table_reference.json was captured from the original 4500-line
if/elif implementation (pre-refactor, git history) by probing every element
that appeared in the chain x every photoelectron line string (plus a
'__default__' sentinel that exercised each element's else branch, and an
unknown-element control). The refactored data-table implementation must
reproduce it exactly.
"""

import json
import pathlib

import pytest

from PhysicsModules.XPS.xps_neo.periodic_table import ElementData, LineParams, get_line_params

REFERENCE = json.loads(
    (
        pathlib.Path(__file__).resolve().parents[1]
        / "component"
        / "references"
        / "periodic_table_reference.json"
    ).read_text()
)

RETURN_NAMES = [
    "BE_lit",
    "is_singlet",
    "so_split",
    "branching_ratio",
    "BE_alt",
    "alt_width",
    "width_range",
    "width",
    "rec_width",
    "Default",
    "peakTypes",
    "is_coster_kronig",
]


def test_reference_is_complete():
    # 95 harvested elements + 1 unknown control, 16 lines + '__default__' each
    assert len(REFERENCE) == 96 * 17


def test_getparams_matches_frozen_reference():
    ed = ElementData(["H"] * 10, ["1s"] * 10)
    mismatches = []
    for key, expected in REFERENCE.items():
        element, line = key.split("|")
        ret = ed.getParams([element] * 10, [line] * 10)
        for name, values in zip(RETURN_NAMES, ret):
            got = values if name == "Default" else values[0]
            if got != expected[name]:
                mismatches.append((key, name, expected[name], got))
    assert not mismatches, mismatches[:20]


def test_getparams_slots_are_independent():
    ed = ElementData(["H"] * 10, ["1s"] * 10)
    elements = ["H", "C", "O", "Ti", "Fe", "Au", "Ta", "Si", "Xx", "N"]
    lines = ["1s", "1s", "2s", "2p", "2p", "4f", "4f", "2p", "1s", "1s"]
    ret = ed.getParams(elements, lines)
    for i in range(10):
        expected = REFERENCE[f"{elements[i]}|{lines[i]}"]
        for name, values in zip(RETURN_NAMES, ret):
            if name == "Default":
                continue
            assert values[i] == expected[name], (i, name)


@pytest.mark.parametrize(
    "element,line,field,value",
    [
        ("Ti", "2p", "is_singlet", False),
        ("Ti", "2p", "so_split", 6.17),
        ("Ti", "2p", "branching_ratio", 0.5),
        ("Ti", "2p", "peak_type", "Double Lorentzian"),
        ("Ti", "2p", "be_lit", 454.1),
        ("Ta", "4f", "is_coster_kronig", True),
        ("C", "1s", "be_lit", 284.5),
        ("Au", "4f", "is_singlet", False),
    ],
)
def test_known_literature_values(element, line, field, value):
    assert getattr(get_line_params(element, line), field) == value


def test_unknown_element_gets_defaults():
    assert get_line_params("Xx", "1s") == LineParams()
    assert get_line_params("Xx", "1s").be_lit == 0.0


def test_unlisted_line_falls_back_to_element_default():
    # Ti has no '5d' branch; the original fell to its else (2p doublet values)
    p = get_line_params("Ti", "5d")
    assert p.be_lit == 454.1
    assert p.is_singlet is False
    assert p.so_split == 6.17
