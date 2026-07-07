"""
Base solver: owns the generation loop and lifecycle hooks. Concrete
solvers (GA, DE, ...) implement step().
"""

from abc import ABC, abstractmethod
from typing import Dict

from Solvers.core.population import Population
from Solvers.core.problem import OptimizationProblem
from Solvers.core.result import SolverResult
from Solvers.core.run_state import RunState


class BaseSolver(ABC):
    name = "base"

    def __init__(self, problem: OptimizationProblem, options: Dict = None, logger=None):
        self.problem = problem
        self.options = dict(options or {})
        self.logger = logger

        self.n_pops = int(self.options.get("nPops", 100))
        self.n_gen = int(self.options.get("nGen", 100))

        self.state = RunState(n_gen=self.n_gen)
        self.population = Population(problem, self.state, n_pops=self.n_pops)
        self.result = SolverResult()
        self._initialized = False

    def log(self, message: str) -> None:
        if self.logger is not None:
            self.logger.print(message)

    def initialize(self) -> None:
        self.population.initialize_populations()
        self._initialized = True

    @abstractmethod
    def step(self) -> None:
        """Advance the population by one generation (in place)."""

    def run(self) -> SolverResult:
        if not self._initialized:
            self.initialize()

        for _ in range(self.n_gen):
            self.state.start_gen()
            self.step()
            self.result.collect(self.state, self.population)
            self.problem.on_generation_end(self.state, self.population)
            self.state.end_gen()

        self.problem.on_run_end(self.state, self.population)
        return self.result

    def __str__(self) -> str:
        return f"{type(self).__name__}(nPops={self.n_pops}, nGen={self.n_gen})"
