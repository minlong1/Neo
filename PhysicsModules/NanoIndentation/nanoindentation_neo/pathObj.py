"""
Fit models for nanoindentation unloading curves.

Ported from nano_indent (pathObj.py). The genome-facing parameter grids are
built in problem.py; these classes remain the model definitions (and keep
the historical object API for analysis/GUI code).
"""

import random

import numpy as np

# Default (low, high, step) grids per Oliver-Pharr parameter
DEFAULT_A_RANGE = (1e-6, 5e-4, 1e-7)
DEFAULT_HF_RANGE = (0.001, 1100, 1)
DEFAULT_M_RANGE = (0.001, 4, 0.01)


def select_opt(def_range, alt_range):
    """
    Used for selecting ranges of using default range

    def_range = if inputs is selected
    alt_range = if inputs is not selected
    """
    if len(def_range) == 0:
        return np.arange(alt_range[0], alt_range[1], alt_range[2])
    return np.arange(def_range[0], def_range[1], def_range[2])


class OliverPharr:
    """
    Power-Law used for nano-indentation

    y = A(x-hf)^m

    """

    n_params = 3
    param_names = ("A", "h_f", "m")

    def __init__(self, path_range, pars_range):
        self.A_range = select_opt(pars_range["A_range"], DEFAULT_A_RANGE)
        self.hf_range = select_opt(pars_range["hf_range"], DEFAULT_HF_RANGE)
        self.m_range = select_opt(pars_range["m_range"], DEFAULT_M_RANGE)

        self.A = np.random.choice(self.A_range)
        self.h_f = np.random.choice(self.hf_range)
        self.m = np.random.choice(self.m_range)

    @staticmethod
    def func(h, A, h_f, m):
        """Vectorized model evaluation used by the fitness function."""
        return A * (h - h_f) ** m

    def get_A(self):
        return self.A

    def get_hf(self):
        return self.h_f

    def get_m(self):
        return self.m

    def get(self):
        return [self.A, self.h_f, self.m]

    def get_func(self, h):
        return self.A * (h - self.h_f) ** self.m

    def set_A(self, A):
        self.A = A

    def set_hf(self, h_f):
        self.h_f = h_f

    def set_m(self, m):
        self.m = m

    def set(self, A, h_f, m):
        self.set_A(A)
        self.set_hf(h_f)
        self.set_m(m)

    def mutate_A(self, chance):
        if random.random() * 100 < chance:
            self.A = np.random.choice(self.A_range)

    def mutate_hf(self, chance):
        if random.random() * 100 < chance:
            self.h_f = np.random.choice(self.hf_range)

    def mutate_m(self, chance):
        if random.random() * 100 < chance:
            self.m = np.random.choice(self.m_range)

    def mutate(self, chance):
        self.mutate_A(chance)
        self.mutate_hf(chance)
        self.mutate_m(chance)

    def verbose(self):
        return [self.get_A(), self.get_hf(), self.get_m()]


MODELS = {
    "oliverpharr": OliverPharr,
}


def get_model(name: str):
    key = name.strip().lower()
    if key not in MODELS:
        raise ValueError(
            f"Unknown fit model {name!r}; supported: {sorted(MODELS)}"
        )
    return MODELS[key]
