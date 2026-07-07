"""
The contract between a physics module and the Solvers package.

A physics module (EXAFS, XPS, NanoIndentation, ...) exposes its fitting task
as an OptimizationProblem: a ParameterSpace describing the genome plus a
fitness function scoring a genome against experimental data. Everything
physics-specific (spectra, models, output files) stays behind this interface.
"""

from abc import ABC, abstractmethod

import numpy as np

from Solvers.core.parameter_space import ParameterSpace


class OptimizationProblem(ABC):
    """Base class every physics module implements to plug into a solver."""

    name: str = "problem"

    def __init__(self, space: ParameterSpace):
        self.space = space

    @abstractmethod
    def fitness(self, genes: np.ndarray) -> float:
        """Score a genome; lower is better."""

    def sample_genes(self) -> np.ndarray:
        """Draw a fresh genome. Override to bias initialization."""
        return self.space.sample()

    # ---- solver lifecycle hooks (all optional) ----

    def on_generation_end(self, state, population) -> None:
        """Called after each generation is evaluated (e.g. write outputs,
        problem-specific parameter sweeps)."""

    def on_run_end(self, state, population) -> None:
        """Called once after the final generation."""
