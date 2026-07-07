"""
Synthetic nanoindentation fixtures: a load-unload displacement curve whose
unloading segment follows a known Oliver-Pharr law, so tests can verify
preprocessing and fitting against ground truth without instrument data.
"""

import numpy as np

# Ground-truth Oliver-Pharr parameters for the synthetic unloading curve
TRUE_A = 2e-4
TRUE_HF = 300.0
TRUE_M = 1.5


def make_curve(n_load=50, n_unload=200, h_max=1000.0):
    """Full load-unload curve: parabolic loading, Oliver-Pharr unloading."""
    h_load = np.linspace(0.0, h_max, n_load)
    y_load = (h_load / h_max) ** 2 * TRUE_A * (h_max - TRUE_HF) ** TRUE_M

    h_unload = np.linspace(h_max, TRUE_HF + 1.0, n_unload)
    y_unload = TRUE_A * (h_unload - TRUE_HF) ** TRUE_M

    h = np.concatenate([h_load, h_unload])
    y = np.concatenate([y_load, y_unload])
    return np.column_stack([h, y])


def write_csv(path, data=None):
    if data is None:
        data = make_curve()
    np.savetxt(path, data, delimiter=",")
    return path
