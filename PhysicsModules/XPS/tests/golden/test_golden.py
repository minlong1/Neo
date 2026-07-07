"""Golden-master matrix as pytest (design.md §3.1).

Run with:  pytest -m golden
Each case shells out to the installed `xps_neo` CLI with a fixed seed and
compares canonicalized outputs bit-exactly against expected/.
"""

import pathlib
import sys

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from run_case import CASES_DIR, EXPECTED_DIR, run_case  # noqa: E402

CASES = sorted(p.stem for p in CASES_DIR.glob("*.ini"))


@pytest.mark.golden
@pytest.mark.parametrize("case", CASES)
def test_golden(case):
    outputs = run_case(case)
    for name, text in outputs.items():
        expected_file = EXPECTED_DIR / case / name
        assert expected_file.exists(), (
            f"no baseline for {case}/{name}; run "
            f"`python tests/golden/run_case.py --capture {case}` and add a "
            f"CHANGELOG.md entry"
        )
        assert text == expected_file.read_text(), (
            f"{case}/{name} diverged from the golden baseline. If this "
            f"change is intended, re-capture and record the reason in "
            f"tests/golden/CHANGELOG.md (design.md §3.1)"
        )
