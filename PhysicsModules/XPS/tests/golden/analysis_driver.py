#!/usr/bin/env python3
"""Headless driver for the GUI post-analysis (gui/xps_analysis2.py).

Replays a completed CLI run directory through the same call sequence the
Tk GUI uses (xps_plot.Analysis_plot.initial_parameters, minus the canvas):
construct -> extract_data -> fitness -> analyze -> get_params, and writes a
canonical, deterministic text dump of every numeric result to out_file.

Usage: python analysis_driver.py <run_dir> <case_ini> <out_file>

Run this in a subprocess: it prepends gui/ to sys.path and imports the
gui-local module names (xps_data, xps_fit, ...), which must not leak into
the caller's process.
"""

import os
import pathlib
import sys

os.environ.setdefault("MPLBACKEND", "Agg")

REPO = pathlib.Path(__file__).resolve().parents[2]

import numpy as np  # noqa: E402

from PhysicsModules.XPS.xps_neo.parser2 import read_input_file  # noqa: E402
from PhysicsModules.XPS.xps_neo.ini_parser import load_config  # noqa: E402


def fmt(value):
    if isinstance(value, np.ndarray):
        value = value.tolist()
    if isinstance(value, (list, tuple)):
        return "[" + ", ".join(fmt(v) for v in value) + "]"
    if isinstance(value, (float, np.floating)):
        return repr(float(value))
    if isinstance(value, (int, np.integer, bool, np.bool_)):
        return (
            repr(int(value))
            if not isinstance(value, (bool, np.bool_))
            else repr(bool(value))
        )
    return repr(str(value))


def main(run_dir, case_ini, out_file):
    config = load_config(read_input_file(case_ini))

    from PhysicsModules.XPS.xps_neo import xps_data  # unified loader (Phase 3c)
    from PhysicsModules.XPS.xps_neo.gui import xps_analysis2  # packaged GUI analysis (Phase 5)

    data_file = os.path.join(run_dir, "data", "DeltaPux1000.txt")
    data_obj = xps_data.xps_data(
        data_file, config["skipLn"], config["x_offset"], config["y_offset"]
    )
    params = {
        "fileName": data_file,
        "peaks": list(config["peak_type"]),
        "data obj": data_obj,
    }

    analysis = xps_analysis2.xps_analysis(run_dir, params, "Voigt")
    analysis.extract_data(False, plot_err=False)
    loss = analysis.fitness(analysis.best_ind)
    analysis.analyze()
    (
        parameters,
        errors,
        errors_bkgns,
        BTou2,
        BTou3,
        peak_areas,
        FWHM_values,
        peak_y_vals,
        totalFit,
        background_fit,
        residual,
        y_raw,
        upper_err_area,
        lower_err_area,
    ) = analysis.get_params()

    lines = []
    lines.append("bestFit: " + fmt(analysis.bestFit))
    lines.append("err: " + fmt(analysis.err))
    lines.append("loss: " + fmt(loss))
    lines.append("parameters: " + fmt(parameters))
    lines.append("errors: " + fmt(errors))
    lines.append("errors_bkgns: " + fmt(errors_bkgns))
    lines.append("BTou2: " + fmt(BTou2))
    lines.append("BTou3: " + fmt(BTou3))
    lines.append("peak_areas: " + fmt(peak_areas))
    lines.append("upper_err_area: " + fmt(upper_err_area))
    lines.append("lower_err_area: " + fmt(lower_err_area))
    lines.append("FWHM_values: " + fmt(FWHM_values))
    lines.append("totalArea: " + fmt(analysis.totalArea))
    lines.append("y_model: " + fmt(totalFit))
    lines.append("background_fit: " + fmt(background_fit))
    for i, comp in enumerate(peak_y_vals):
        lines.append(f"peak_component_{i}: " + fmt(comp))
    pathlib.Path(out_file).write_text("\n".join(lines) + "\n")
    print(f"analysis dump written: {out_file}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
