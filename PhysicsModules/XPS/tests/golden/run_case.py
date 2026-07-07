#!/usr/bin/env python3
"""Golden-master runner for XPS_Neo (design.md §3.1).

Runs the installed `xps_neo` CLI on each case INI in an isolated temp
directory with a fixed seed, canonicalizes the outputs, and either
captures them as the expected baseline or checks them against it.

Usage:
    python tests/golden/run_case.py --check  [case ...]   # default: all
    python tests/golden/run_case.py --capture [case ...]  # (re)baseline

Capturing over existing baselines is a deliberate act: record the reason
in tests/golden/CHANGELOG.md (see design.md §3.1 "Regeneration").
"""

import argparse
import difflib
import pathlib
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from canonicalize import canonicalize_run  # noqa: E402
from tol_compare import compare as tol_compare  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parents[1]
CASES_DIR = HERE / "cases"
EXPECTED_DIR = HERE / "expected"
DATA_SRC = REPO / "datafile" / "DeltaPux1000.txt"
SEED = "42"
TIMEOUT_S = 600


def run_case(case):
    """Run one case in a temp dir; return canonical outputs dict."""
    ini = CASES_DIR / f"{case}.ini"
    if not ini.exists():
        raise FileNotFoundError(ini)
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
        # The CLI exits 0 even on internal errors (bug register #6),
        # so also require that generation lines were actually written.
        gen_file = tmp / "fit_out.txt"
        ran = gen_file.exists() and len(gen_file.read_text().splitlines()) > 1
        if proc.returncode != 0 or not ran:
            raise RuntimeError(
                f"case {case!r} did not complete "
                f"(rc={proc.returncode}, output lines missing)\n"
                f"--- stdout tail ---\n{proc.stdout[-2000:]}\n"
                f"--- stderr tail ---\n{proc.stderr[-2000:]}"
            )
        return canonicalize_run(tmp)


def capture(case):
    outputs = run_case(case)
    dest = EXPECTED_DIR / case
    dest.mkdir(parents=True, exist_ok=True)
    for name, text in outputs.items():
        (dest / name).write_text(text)
    print(f"CAPTURED {case}")


def check(case, rtol=None):
    outputs = run_case(case)
    dest = EXPECTED_DIR / case
    failures = []
    for name, text in outputs.items():
        expected_file = dest / name
        if not expected_file.exists():
            failures.append(f"{name}: no expected file (run --capture first)")
            continue
        expected = expected_file.read_text()
        if rtol is not None and text != expected:
            ok, worst, msg = tol_compare(expected, text, rtol=rtol)
            if ok:
                print(f"  {case}/{name}: {msg}")
                continue
            failures.append(f"{name}: {msg}")
            continue
        if text != expected:
            diff = "\n".join(
                difflib.unified_diff(
                    expected.splitlines(),
                    text.splitlines(),
                    f"expected/{case}/{name}",
                    f"actual/{case}/{name}",
                    lineterm="",
                    n=1,
                )
            )
            failures.append(f"{name} differs:\n{diff[:3000]}")
    if failures:
        print(f"FAIL {case}")
        for f in failures:
            print(f"  {f}")
        return False
    print(f"PASS {case}")
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--capture", action="store_true")
    mode.add_argument("--check", action="store_true")
    ap.add_argument(
        "--rtol",
        type=float,
        default=None,
        help="tolerance-aware check (design.md Phase 4 order-changing "
        "class); prints the worst deviation per file",
    )
    ap.add_argument(
        "cases", nargs="*", help="case names (default: every .ini under cases/)"
    )
    args = ap.parse_args()

    cases = args.cases or sorted(p.stem for p in CASES_DIR.glob("*.ini"))
    ok = True
    for case in cases:
        try:
            if args.capture:
                capture(case)
            else:
                ok = check(case, rtol=args.rtol) and ok
        except Exception as exc:  # surface which case broke, keep going
            print(f"ERROR {case}: {exc}")
            ok = False
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
