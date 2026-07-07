"""
Differential Evolution solver.

Carried over as a stub from EXAFS neoSolver.py (NeoSolver_DE was
unimplemented there too). Fill in step() with mutation/crossover/selection
on gene vectors: for each individual x, build a trial vector from
a + F * (b - c) and keep the better of trial vs x.
"""

from Solvers.core.base_solver import BaseSolver


class DESolver(BaseSolver):
    name = "DE"

    def step(self):
        # TODO: implement differential evolution
        raise NotImplementedError("DE solver is not implemented yet")
