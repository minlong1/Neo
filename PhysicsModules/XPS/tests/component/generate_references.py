#!/usr/bin/env python3
"""(Re)generate .npy reference curves for the component tests.

Regenerating is a rebaselining act — record the reason in
tests/golden/CHANGELOG.md, same policy as the golden matrix.

Usage: python tests/component/generate_references.py
"""

import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import shapes  # noqa: E402
from conftest import load_xps_module  # noqa: E402

REF_DIR = HERE / "references"


def main():
    REF_DIR.mkdir(exist_ok=True)
    xf = load_xps_module("PhysicsModules.XPS.xps_neo.xps_fit")
    for name in shapes.PEAK_VARIANTS:
        y = shapes.eval_peak(xf, name)
        np.save(REF_DIR / f"peak_{name}.npy", y)
        print(
            f"peak_{name}: max={y.max():.6g} at BE="
            f"{shapes.X_GRID[np.argmax(y)]:.2f}, nan={np.isnan(y).any()}"
        )
    for name in shapes.BKGN_VARIANTS:
        y = shapes.eval_background(xf, name)
        np.save(REF_DIR / f"bkgn_{name}.npy", y)
        print(
            f"bkgn_{name}: range=[{y.min():.6g}, {y.max():.6g}], "
            f"nan={np.isnan(y).any()}"
        )


if __name__ == "__main__":
    main()
