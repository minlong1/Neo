"""
Physics-agnostic population-based solvers for the NEO framework.

A physics module plugs in by implementing Solvers.core.OptimizationProblem
(a ParameterSpace plus a fitness function), then running any registered
solver against it:

    from Solvers import get_solver
    solver = get_solver("GA")(problem, options={"nGen": 100, "nPops": 100})
    result = solver.run()

This package must stay free of physics dependencies (no larch, no module-
specific imports) — numpy only.
"""

from Solvers.core import (
    BaseSolver,
    ContinuousGeneRange,
    GeneRange,
    Individual,
    OptimizationProblem,
    ParameterSpace,
    Population,
    RunState,
    SolverResult,
)
from Solvers.ga import GASolver, GARechenbergSolver
from Solvers.de import DESolver
from Solvers.pso import PSOSolver
from Solvers.cmaes import CMAESSolver
from Solvers.demcmc import DEMCMCSolver

# Numeric IDs preserve the historical EXAFS solOpt/solver_type values.
SOLVER_REGISTRY = {
    "GA": GASolver,
    "GA_RECHENBERG": GARechenbergSolver,
    "DE": DESolver,
    "PSO": PSOSolver,
    "CMA_ES": CMAESSolver,
    "DE_MCMC": DEMCMCSolver,
    0: GASolver,
    1: GARechenbergSolver,
    2: DESolver,
}


def get_solver(name_or_id):
    """Look up a solver class by name (case-insensitive) or numeric option ID."""
    key = name_or_id.upper() if isinstance(name_or_id, str) else name_or_id
    if key not in SOLVER_REGISTRY:
        raise ValueError(
            f"Unknown solver {name_or_id!r}; available: "
            f"{sorted(k for k in SOLVER_REGISTRY if isinstance(k, str))}"
        )
    return SOLVER_REGISTRY[key]


__all__ = [
    "BaseSolver",
    "ContinuousGeneRange",
    "GeneRange",
    "Individual",
    "OptimizationProblem",
    "ParameterSpace",
    "Population",
    "RunState",
    "SolverResult",
    "GASolver",
    "GARechenbergSolver",
    "DESolver",
    "PSOSolver",
    "CMAESSolver",
    "DEMCMCSolver",
    "SOLVER_REGISTRY",
    "get_solver",
]
