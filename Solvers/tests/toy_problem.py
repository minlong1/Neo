"""
Analytic toy problem used by the Solvers test suite: minimize
sum((x - target)^2) over a discrete grid. Global optimum is the grid
point closest to `target` in every dimension.
"""

import numpy as np

from Solvers.core import GeneRange, OptimizationProblem, ParameterSpace


class QuadraticProblem(OptimizationProblem):
    name = "quadratic"

    def __init__(self, target=None, n_genes: int = 4):
        if target is None:
            target = np.linspace(-0.5, 0.5, n_genes)
        self.target = np.asarray(target, dtype=float)
        space = ParameterSpace(
            [
                GeneRange(np.arange(-1.0, 1.001, 0.01), name=f"x{i}")
                for i in range(len(self.target))
            ]
        )
        super().__init__(space)
        self.generation_end_calls = 0
        self.run_end_calls = 0

    def fitness(self, genes: np.ndarray) -> float:
        return float(np.sum((np.asarray(genes) - self.target) ** 2))

    def on_generation_end(self, state, population) -> None:
        self.generation_end_calls += 1

    def on_run_end(self, state, population) -> None:
        self.run_end_calls += 1
