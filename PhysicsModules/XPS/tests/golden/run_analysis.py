#!/usr/bin/env python3
"""Golden-master runner for the GUI post-analysis (design.md Phase 3).

For each case: runs the CLI with the fixed seed into a temp dir (exactly
like run_case.py), then replays that run directory through
gui/xps_analysis2 via analysis_driver.py in a subprocess, and captures or
checks the canonical analysis dump under expected_analysis/<case>.txt.

Usage:
    python tests/golden/run_analysis.py --check   [case ...]
    python tests/golden/run_analysis.py --capture [case ...]

Recapturing is a rebaselining act: record the reason in CHANGELOG.md.
"""

import argparse
import difflib
import pathlib
import shutil
import subprocess
import sys
import tempfile

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from run_case import CASES_DIR, DATA_SRC, SEED, TIMEOUT_S  # noqa: E402
from tol_compare import compare as tol_compare  # noqa: E402

EXPECTED_DIR = HERE / "expected_analysis"
DRIVER = HERE / "analysis_driver.py"

# Analysis-golden subset: one case per analysis code path that works today.
# Excluded: Doniach_shirley and voigt_coster_kronig (the analysis error-splitting
# for these is flagged buggy in gui code comments), voigt_de_mode
# (identical analysis path to voigt_baseline_svsc), lorentzian_baseline
# (same path as gaussian_linear apart from the peak formula, which the
# component tests already pin).

CASES = [
    "voigt_baseline_svsc",
    "voigt_shirley",
    "gaussian_linear",
    "double_lorentzian_doublet",
    "two_peak_limited_correlated",
]


def run_analysis(case):
    ini = CASES_DIR / f"{case}.ini"
    with tempfile.TemporaryDirectory() as tmp:
        tmp = pathlib.Path(tmp)
        shutil.copy(ini, tmp / "case.ini")
        (tmp / "data").mkdir()
        shutil.copy(DATA_SRC, tmp / "data" / DATA_SRC.name)
        proc = subprocess.run(
            ["xps_neo", "-i", "case.ini", "--seed", SEED],
            cwd=tmp,
            capture_output=True,
            text=True,
            timeout=TIMEOUT_S,
        )
        gen_file = tmp / "fit_out.txt"
        if (
            proc.returncode != 0
            or not gen_file.exists()
            or len(gen_file.read_text().splitlines()) <= 1
        ):
            raise RuntimeError(
                f"CLI run for {case!r} failed:\n"
                f"{proc.stdout[-1500:]}\n{proc.stderr[-1500:]}"
            )
        out = tmp / "analysis_dump.txt"
        proc = subprocess.run(
            [sys.executable, str(DRIVER), str(tmp), str(tmp / "case.ini"), str(out)],
            capture_output=True,
            text=True,
            timeout=TIMEOUT_S,
        )
        if proc.returncode != 0 or not out.exists():
            raise RuntimeError(
                f"analysis driver for {case!r} failed:\n"
                f"{proc.stdout[-1500:]}\n{proc.stderr[-1500:]}"
            )
        return out.read_text()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--capture", action="store_true")
    mode.add_argument("--check", action="store_true")
    ap.add_argument("--rtol", type=float, default=None)
    ap.add_argument("cases", nargs="*")
    args = ap.parse_args()

    ok = True
    for case in args.cases or CASES:
        try:
            text = run_analysis(case)
            expected_file = EXPECTED_DIR / f"{case}.txt"
            if args.capture:
                EXPECTED_DIR.mkdir(exist_ok=True)
                expected_file.write_text(text)
                print(f"CAPTURED {case}")
            elif args.rtol is not None and text != expected_file.read_text():
                ok2, worst, msg = tol_compare(
                    expected_file.read_text(), text, rtol=args.rtol
                )
                if ok2:
                    print(f"PASS {case} ({msg})")
                else:
                    print(f"FAIL {case}: {msg}")
                    ok = False
            else:
                expected = expected_file.read_text()
                if text != expected:
                    diff = "\n".join(
                        difflib.unified_diff(
                            expected.splitlines(),
                            text.splitlines(),
                            "expected",
                            "actual",
                            lineterm="",
                            n=0,
                        )
                    )
                    print(f"FAIL {case}\n{diff[:2500]}")
                    ok = False
                else:
                    print(f"PASS {case}")
        except Exception as exc:
            print(f"ERROR {case}: {exc}")
            ok = False
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
