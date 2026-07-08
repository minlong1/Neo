"""
Generic solver result: per-generation history plus the best individual.
Ported from the EXAFS NeoResult minus physics-specific statistics.
"""

import pickle
from typing import List


class SolverResult:
    """Post-run summary: the best genome found and its per-generation
    trajectory (`historyBest` is the running global best; `genBest` is that
    generation's own best, which can be worse)."""

    def __init__(self):
        self.best_individual = None
        self.best_value: float = None
        self.historyBest: List[float] = []
        self.genBest: List[float] = []

    def collect(self, state, population) -> None:
        self.best_individual = state.globBestInd
        self.best_value = state.globBestVal
        self.historyBest.append(state.globBestVal)
        self.genBest.append(state.currBestVal)

    def save(self, filename: str) -> None:
        with open(filename, "wb") as file:
            pickle.dump(self, file)

    @classmethod
    def load(cls, filename: str) -> "SolverResult":
        with open(filename, "rb") as file:
            return pickle.load(file)
        return None

    def __str__(self) -> str:
        if self.best_individual is None:
            return "Best Individual: None"
        return f"Best Individual: {self.best_individual}, value: {self.best_value}"
