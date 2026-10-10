"""
Raman spectrum loading and preprocessing.

Raman analogue of NanoIndentation's ``nano_neo_data.py`` (and of larch's
``read_ascii`` + ``autobk`` in the EXAFS module): read a two-column
(Raman shift, intensity) ascii spectrum, resample it onto a uniform model
grid, optionally subtract a polynomial baseline, and normalise to unit peak
inside the fit window.

Grid convention (unchanged from raman_neo): ``xgrid = arange(big + 1) * dx``
spans 0..xmax, and ``interval`` holds the indices of the fit window
``[xmin, xmax]``.
"""

import numpy as np

from PhysicsModules.Raman.raman_neo.bandshape import polynomial_baseline


def read_spectrum(path):
    """Load a two-column (wavenumber, intensity) ascii spectrum, sorted by x.

    Comma- or whitespace-delimited; lines starting with ``#`` are comments.
    """
    with open(path) as fh:
        sample = "".join(line for line in fh if line.strip() and not line.lstrip().startswith("#"))
    delimiter = "," if "," in sample else None
    data = np.atleast_2d(np.genfromtxt(str(path), comments="#", delimiter=delimiter))
    if data.shape[1] < 2:
        raise ValueError(f"Data file {path} must have at least two columns (wavenumber, intensity)")
    x = data[:, 0].astype(float)
    y = data[:, 1].astype(float)
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    order = np.argsort(x)
    return x[order], y[order]


def build_grid(xmin, xmax, dx):
    """Model grid 0..xmax with spacing dx, plus the fit-window indices."""
    if xmax <= xmin:
        raise ValueError(f"xmax ({xmax}) must be greater than xmin ({xmin})")
    if dx <= 0:
        raise ValueError(f"deltax must be positive, got {dx}")
    small = int(xmin / dx)
    big = int(xmax / dx)
    xgrid = np.arange(big + 1) * dx
    interval = np.arange(small, big + 1)
    return xgrid, interval


class RamanSpectrum:
    """An experimental spectrum resampled onto the model grid.

    Attributes: ``raw_x``/``raw_y`` (file as read), ``xgrid``/``interval``
    (see module docstring), ``exp`` (processed intensity on ``xgrid``),
    ``baseline`` (subtracted baseline on ``xgrid``), and ``x``/``y`` (the
    fit window: ``xgrid[interval]``/``exp[interval]``).
    """

    def __init__(
        self,
        data_file,
        xmin: float = 100.0,
        xmax: float = 1800.0,
        dx: float = 1.0,
        baseline_order: int = 0,
        normalize: bool = True,
    ):
        self.data_file = data_file
        self.xmin = float(xmin)
        self.xmax = float(xmax)
        self.dx = float(dx)
        self.baseline_order = int(baseline_order)
        self.normalize = bool(normalize)

        self.raw_x, self.raw_y = read_spectrum(data_file)
        self.xgrid, self.interval = build_grid(self.xmin, self.xmax, self.dx)
        self.exp, self.baseline = preprocess(
            self.raw_x,
            self.raw_y,
            self.xgrid,
            baseline_order=self.baseline_order,
            normalize=self.normalize,
            interval=self.interval,
        )

    @property
    def x(self) -> np.ndarray:
        return self.xgrid[self.interval]

    @property
    def y(self) -> np.ndarray:
        return self.exp[self.interval]

    def get_raw_data(self) -> np.ndarray:
        return np.column_stack([self.raw_x, self.raw_y])

    def get_fit_data(self) -> np.ndarray:
        return np.column_stack([self.x, self.y])


def preprocess(raw_x, raw_y, xgrid, baseline_order=0, normalize=True, interval=None):
    """Resample onto ``xgrid``, subtract a baseline, normalise to unit peak.

    ``raw_x`` must be sorted ascending (``read_spectrum`` does this).

    Returns ``(exp, baseline)``, both on ``xgrid``. Normalisation (and the
    peak it divides by) is taken within ``interval`` when given.
    """
    exp = np.interp(xgrid, raw_x, raw_y, left=0.0, right=0.0)

    baseline = np.zeros_like(exp)
    if baseline_order and baseline_order > 0:
        # Fit only where the spectrum has data: grid points outside it are
        # zero padding, which would drag the lower-envelope fit down.
        covered = (xgrid >= raw_x[0]) & (xgrid <= raw_x[-1])
        if np.count_nonzero(covered) > baseline_order:
            baseline[covered] = polynomial_baseline(xgrid[covered], exp[covered], order=baseline_order)
            exp = exp - baseline

    if normalize:
        idx = interval if interval is not None else slice(None)
        window = exp[idx]
        peak = np.max(window) if window.size and np.max(window) > 0 else 1.0
        exp = exp / peak
        baseline = baseline / peak

    return exp, baseline
