"""
Stagnation detection and hyper-mutation (Forbes & Long, DE-BNN, Section 3).

The "residual" is the best-fitness change from one generation to the next
(negative = improvement). Stagnation is declared when the t-generation
running average of the residual is ~0 — the search has stopped improving.
While stagnant, F, CR, and the mutation operator are resampled each
generation from the paper's hyper-mutation ranges (Eq. 11, 12, 15) instead
of using the fixed defaults, to diversify the search; once a strictly
better candidate is found (a negative residual breaks the streak), the
schedule reverts to the defaults.
"""

from typing import List, Optional, Tuple

import numpy as np

# Eq. 11: F_i ~ U{0.1, 0.9} at 0.01 intervals.
F_CHOICES = np.round(np.arange(0.10, 0.901, 0.01), 2)
# Eq. 15: CR_i ~ U{0.5, 0.9}. (Section 3.4's prose describes a 0.10-0.9-by-0.1
# range, apparently reusing Section 3.1's F wording verbatim; Eq. 15 itself
# is explicit, so we use its literal {0.5, 0.6, 0.7, 0.8, 0.9}.)
CR_CHOICES = np.round(np.arange(0.5, 0.901, 0.1), 2)
# Eq. 12: mutation operator drawn from {DE/best/j, DE/rand/j}; this
# implementation supports j in {1, 2} (see mutation.py).
OPERATOR_CHOICES = ["best/1", "best/2", "rand/1", "rand/2"]


class StagnationTracker:
    """Tracks a running-average residual over `window` generations and
    flags stagnation when it's within `tol` of zero (see module docstring
    for what "residual" means here)."""

    def __init__(self, window: int = 10, tol: float = 1e-9):
        self.window = window
        self.tol = tol
        self._residuals: List[float] = []
        self._prev_best: Optional[float] = None
        self.iteration = 0  # generations spent in the current stagnation streak

    def update(self, current_best: float) -> bool:
        """Record this generation's best fitness; return whether stagnant."""
        if self._prev_best is not None:
            residual = current_best - self._prev_best
            self._residuals.append(residual)
            if len(self._residuals) > self.window:
                self._residuals.pop(0)
        self._prev_best = current_best

        stagnant = self.is_stagnant()
        self.iteration = self.iteration + 1 if stagnant else 0
        return stagnant

    def is_stagnant(self) -> bool:
        if len(self._residuals) < self.window:
            return False
        return abs(float(np.mean(self._residuals))) < self.tol


def sample_hyperparameters() -> Tuple[float, float, str]:
    """Draw a fresh (F, CR, mutation_operator) triple for hyper-mutation."""
    F = float(np.random.choice(F_CHOICES))
    CR = float(np.random.choice(CR_CHOICES))
    operator = str(np.random.choice(OPERATOR_CHOICES))
    return F, CR, operator
