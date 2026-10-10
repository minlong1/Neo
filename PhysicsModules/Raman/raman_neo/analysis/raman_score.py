"""
Standalone scoring helpers, independent of the solvers.

Raman analogue of the EXAFS module's ``analysis/larch_score.py``: given band
definitions and parameters, synthesise a model spectrum and score it against
an experimental one. Useful for scoring a fixed candidate model without
running a fit. Spectrum preprocessing is shared with RamanProblem
(``raman_data.preprocess``), so scores here match the fit's fitness.
"""

import numpy as np

from PhysicsModules.Raman.raman_neo.bandshape import Band, read_band_file
from PhysicsModules.Raman.raman_neo.raman_data import build_grid, preprocess, read_spectrum

__all__ = [
    "build_grid",
    "load_experimental",
    "synthesize",
    "score",
    "r_factor",
    "bands_from_file",
]


def load_experimental(path, xgrid, baseline_order=0, normalize=True, interval=None):
    """Load, resample, baseline-correct and (optionally) normalise a spectrum."""
    x, y = read_spectrum(path)
    exp, _ = preprocess(x, y, xgrid, baseline_order=baseline_order, normalize=normalize, interval=interval)
    return exp


def synthesize(bands, params, xgrid):
    """Total spectrum from Band objects and an (nband, 4) param array whose
    rows are ``[amp, shift, width, center_shift]``."""
    params = np.asarray(params, dtype=float)
    y_total = np.zeros_like(xgrid, dtype=float)
    for band, row in zip(bands, params):
        band.set_params(*row)
        y_total += band.band2chi(xgrid)
    return y_total


def score(y_model, y_exp, xgrid=None, interval=None, weight=0.0):
    """Weighted sum-of-squares misfit (the fit's fitness, computed standalone)."""
    if interval is None:
        interval = np.arange(len(y_model))
    diff = y_model[interval] - y_exp[interval]
    if weight and weight != 0.0 and xgrid is not None:
        diff = diff * (xgrid[interval] ** weight)
    return float(np.sum(diff**2))


def r_factor(y_model, y_exp, interval=None):
    """EXAFS-style R-factor: sum |model-exp|^2 / sum |exp|^2."""
    if interval is None:
        interval = np.arange(len(y_model))
    num = np.sum((y_model[interval] - y_exp[interval]) ** 2)
    den = np.sum(y_exp[interval] ** 2)
    return float(num / den) if den > 0 else float("inf")


def bands_from_file(band_file, band_list=None):
    """Instantiate Band objects from a band-definition file (1-based rows)."""
    defs = read_band_file(band_file)
    if band_list is None:
        band_list = list(range(1, len(defs) + 1))
    bands = []
    for band_id in band_list:
        d = defs[int(band_id) - 1]
        bands.append(Band(f"Band{band_id}", d["center"], d["shape"], d["eta"], d["base_width"]))
    return bands
