"""
Synthetic Raman fixtures: a noise-free spectrum built from known bands, so
tests can verify preprocessing, fitness, and fitting against ground truth
without instrument data.
"""

import numpy as np

from PhysicsModules.Raman.raman_neo.bandshape import lorentzian, pseudo_voigt

# Ground-truth bands: (center, amplitude, width, shape, eta)
TRUE_BANDS = [
    (480.0, 0.45, 9.0, "lorentzian", 1.0),
    (720.0, 0.30, 12.0, "lorentzian", 1.0),
    (1350.0, 0.80, 28.0, "voigt", 0.7),
    (1580.0, 1.00, 18.0, "lorentzian", 1.0),
]
TRUE_SHIFT = 0.0
# Decoy rows in the band file with no real intensity behind them
DECOY_BANDS = [
    (1000.0, "lorentzian", 1.0, 15.0),
    (1700.0, "gaussian", 0.0, 20.0),
]


def make_spectrum(x=None, shift=TRUE_SHIFT, baseline=None, noise=0.0, seed=0):
    """Two-column (x, intensity) array of the ground-truth bands."""
    if x is None:
        x = np.arange(100.0, 1800.0 + 1.0, 1.0)
    y = np.zeros_like(x)
    for center, amp, width, shape, eta in TRUE_BANDS:
        c = center + shift
        if shape == "voigt":
            y += amp * pseudo_voigt(x, c, width, eta)
        else:
            y += amp * lorentzian(x, c, width)
    if baseline is not None:
        y += baseline(x)
    if noise:
        y += np.random.default_rng(seed).normal(0.0, noise, size=x.shape)
    return np.column_stack([x, y])


def true_genes(nband=4, shift=TRUE_SHIFT):
    """Flat genome [shift, (amp, width, cshift) x nband] of the truth; extra
    bands beyond the real four get zero amplitude."""
    genes = [shift]
    for i in range(nband):
        if i < len(TRUE_BANDS):
            _, amp, width, _, _ = TRUE_BANDS[i]
            genes += [amp, width, 0.0]
        else:
            genes += [0.0, 10.0, 0.0]
    return np.array(genes)


def write_spectrum(path, data=None, delimiter=" "):
    if data is None:
        data = make_spectrum()
    np.savetxt(path, data, delimiter=delimiter, header="wavenumber intensity")
    return path


def write_band_file(path, decoys=True):
    lines = ["# center, shape, eta, base_width"]
    for center, _, width, shape, eta in TRUE_BANDS:
        lines.append(f"{center}, {shape}, {eta}, {width}")
    if decoys:
        for center, shape, eta, width in DECOY_BANDS:
            lines.append(f"{center}, {shape}, {eta}, {width}")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path
