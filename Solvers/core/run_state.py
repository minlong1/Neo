"""
Mutable run state shared between solver, population, operators, and
problem hooks: generation counter, timing, and best-fit tracking.
"""

import time

import numpy as np


class RunState:
    def __init__(self, n_gen: int = 100):
        self.nGen = n_gen
        self.currGen = 0

        # Best-fit tracking
        self.currBestInd = None
        self.currBestVal = np.inf
        self.globBestInd = None
        self.globBestVal = np.inf
        self.bestDiff = np.inf

        # Adaptive-mutation bookkeeping (Rechenberg)
        self.diffCounter = 0

        # Timing
        self.tt = 0.0
        self.currGen_st = 0.0
        self.currGen_tt = 0.0

    @property
    def second_half(self) -> bool:
        return self.currGen >= self.nGen // 2

    def start_gen(self) -> None:
        self.currGen += 1
        self.currGen_st = time.time()

    def end_gen(self) -> None:
        self.currGen_tt = time.time() - self.currGen_st
        self.tt += self.currGen_tt

    def update_best(self, best_individual, best_value: float, min_delta: float = 0.01) -> None:
        """Record the generation's best and refresh the global best.

        `bestDiff` is the improvement over the global best, zeroed when the
        change is below `min_delta` (used by Metropolis/Rechenberg schedules).
        """
        self.currBestInd = best_individual
        self.currBestVal = best_value

        delta = self.globBestVal - self.currBestVal
        self.bestDiff = delta if np.abs(delta) > min_delta else 0

        if self.currBestVal < self.globBestVal:
            self.globBestInd = best_individual
            self.globBestVal = best_value
