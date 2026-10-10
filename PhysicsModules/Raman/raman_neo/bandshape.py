"""
Raman band line-shape synthesis.

The Raman analogue of larch's ``feffdat`` in the EXAFS module: it turns fit
parameters into a synthetic spectrum. Instead of photoelectron scattering
paths, each *band* is an analytic line shape (Lorentzian, Gaussian, or
pseudo-Voigt) centred on a nominal Raman shift. Pure numpy.

A band has four fit parameters, mirroring the four EXAFS path parameters:

    ======================  ======================  ==============================
    EXAFS path parameter    Raman band parameter    meaning
    ======================  ======================  ==============================
    S02                     amp                     band amplitude (peak height)
    E0                      shift (global)          global wavenumber offset
    sigma2                  width                   band width (HWHM/sigma, cm^-1)
    deltaR                  center_shift            per-band centre offset (cm^-1)
    ======================  ======================  ==============================

The absolute centre of a band is ``nominal_center + center_shift + shift``.
``shift`` (the E0 analogue) is shared across every band in a genome and acts
like an x-axis calibration offset.

Ported from the standalone ``raman_neo`` package (bandshape.py).
"""

import numpy as np

# numpy>=2 renamed ``trapz`` to ``trapezoid``; support both.
_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))


# --------------------------------------------------------------------------- #
# Line shapes. ``width`` is the Lorentzian HWHM / Gaussian sigma in cm^-1.
# Each shape has unit peak height so ``amp`` is the physical peak intensity.
# --------------------------------------------------------------------------- #
def lorentzian(x, center, width):
    width = max(float(width), 1e-9)
    return width**2 / ((x - center) ** 2 + width**2)


def gaussian(x, center, width):
    width = max(float(width), 1e-9)
    return np.exp(-((x - center) ** 2) / (2.0 * width**2))


def pseudo_voigt(x, center, width, eta=0.5):
    """Linear combination of a Lorentzian and a Gaussian.

    ``eta`` in [0, 1] is the Lorentzian fraction (eta=1 -> pure Lorentzian).
    """
    eta = float(np.clip(eta, 0.0, 1.0))
    return eta * lorentzian(x, center, width) + (1.0 - eta) * gaussian(x, center, width)


SHAPES = {
    "lorentzian": lorentzian,
    "gaussian": gaussian,
    "voigt": pseudo_voigt,
    "pseudo_voigt": pseudo_voigt,
}


def evaluate_shape(shape, x, center, width, eta=0.5):
    shape = str(shape).lower()
    func = SHAPES.get(shape)
    if func is None:
        raise ValueError(f"Unknown line shape {shape!r}; choose from {list(SHAPES)}")
    if func is pseudo_voigt:
        return pseudo_voigt(x, center, width, eta)
    return func(x, center, width)


# --------------------------------------------------------------------------- #
# Band object: the synthesis analogue of a larch ``feffpath``.
# --------------------------------------------------------------------------- #
class Band:
    """A single spectral band: its fixed description (nominal centre, shape,
    Lorentzian fraction) plus the current fit parameters, set before each
    :meth:`band2chi` call (mirrors a larch ``feffpath``)."""

    def __init__(self, name, center, shape="lorentzian", eta=0.5, base_width=10.0):
        if str(shape).lower() not in SHAPES:
            raise ValueError(f"Unknown line shape {shape!r}; choose from {list(SHAPES)}")
        self.name = name
        self.center = float(center)  # nominal centre (cm^-1)
        self.shape = str(shape).lower()
        self.eta = float(eta)
        self.base_width = float(base_width)  # reference width (informational)

        # Fit parameters
        self.amp = 0.0  # <- S02
        self.shift = 0.0  # <- E0 (global, shared)
        self.width = float(base_width)  # <- sigma2
        self.center_shift = 0.0  # <- deltaR

        # Result cache
        self.x = None
        self.chi = None

    def set_params(self, amp, shift, width, center_shift):
        self.amp = float(amp)
        self.shift = float(shift)
        self.width = float(width)
        self.center_shift = float(center_shift)

    def absolute_center(self):
        return self.center + self.center_shift + self.shift

    def evaluate(self, x, amp, shift, width, center_shift):
        """Stateless evaluation (used by the fitness hot path)."""
        center = self.center + center_shift + shift
        return amp * evaluate_shape(self.shape, x, center, width, self.eta)

    def band2chi(self, x):
        """Evaluate this band on grid ``x``. Analogue of ``feffdat.path2chi``."""
        x = np.asarray(x, dtype=float)
        y = self.evaluate(x, self.amp, self.shift, self.width, self.center_shift)
        self.x = x
        self.chi = y
        return y

    def area(self, x):
        """Integrated intensity of the band on grid ``x``."""
        if self.chi is None or self.x is None or not np.array_equal(self.x, x):
            self.band2chi(x)
        return float(_trapz(np.abs(self.chi), x))

    def __str__(self):
        return (
            f"Band {self.name}: center={np.round(self.center, 2)}, "
            f"amp={np.round(self.amp, 3)}, shift={np.round(self.shift, 2)}, "
            f"width={np.round(self.width, 3)}, center_shift={np.round(self.center_shift, 3)}"
        )


# --------------------------------------------------------------------------- #
# Band definition file reader (analogue of the FEFF path folder).
# --------------------------------------------------------------------------- #
def read_band_file(filepath):
    """Read a band-definition file.

    Format: one band per line, comma separated::

        # center, shape, eta, base_width
        520.0, lorentzian, 1.0, 8.0
        1350.0, voigt, 0.5, 25.0

    Only the centre is required; ``shape`` defaults to ``lorentzian``,
    ``eta`` to 0.5 and ``base_width`` to 10.0. Lines starting with ``#`` and
    blank lines are ignored.

    Returns a list of dicts with keys ``center``, ``shape``, ``eta``,
    ``base_width``.
    """
    bands = []
    with open(filepath) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.split(",")]
            center = float(parts[0])
            shape = parts[1].lower() if len(parts) > 1 and parts[1] else "lorentzian"
            eta = float(parts[2]) if len(parts) > 2 and parts[2] else 0.5
            base_width = float(parts[3]) if len(parts) > 3 and parts[3] else 10.0
            bands.append({"center": center, "shape": shape, "eta": eta, "base_width": base_width})
    if not bands:
        raise ValueError(f"No bands found in band file {filepath!r}")
    return bands


# --------------------------------------------------------------------------- #
# Baseline removal (analogue of ``autobk`` background subtraction).
# --------------------------------------------------------------------------- #
def polynomial_baseline(x, y, order=3, niter=40, tol=1e-3):
    """Estimate a smooth baseline via iterative polynomial fitting.

    At each iteration the polynomial fit is clipped to lie below the data
    (points above the fit are pulled down to the fit), which converges to a
    baseline that follows the lower envelope of the spectrum. A light,
    dependency-free stand-in for asymmetric-least-squares baselines.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if order <= 0:
        return np.zeros_like(y)

    # Normalise x for numerical stability of polyfit.
    xs = (x - x.mean()) / (x.std() if x.std() > 0 else 1.0)
    z = y.copy()
    baseline = np.zeros_like(y)
    for _ in range(niter):
        coeffs = np.polyfit(xs, z, order)
        baseline = np.polyval(coeffs, xs)
        new_z = np.minimum(z, baseline)
        if np.linalg.norm(new_z - z) / (np.linalg.norm(z) + 1e-12) < tol:
            z = new_z
            break
        z = new_z
    return baseline
