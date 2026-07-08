"""
Parameter space abstraction for population-based solvers.

A genome is a flat vector of genes. Each gene draws its values from a
GeneRange: a discrete, ordered set of allowed values (mirroring how the
physics modules define fitting-parameter grids, e.g. EXAFS path ranges).

ContinuousGeneRange is a sibling for genes that aren't naturally a
discrete grid — e.g. neural-network weights/biases (Solvers.demcmc), which
are drawn from a continuous initialization distribution (He/Xavier-style)
rather than enumerated candidate values.
"""

from typing import List, Optional, Sequence

import numpy as np


class GeneRange:
    """Discrete set of allowed values for a single gene."""

    def __init__(self, values: Sequence[float], name: str = ""):
        values = np.asarray(values, dtype=float)
        if values.size == 0:
            raise ValueError("GeneRange requires at least one allowed value")
        self.values = values
        self.name = name

    @classmethod
    def from_bounds(cls, low: float, high: float, step: float, name: str = "") -> "GeneRange":
        return cls(np.arange(low, high, step), name=name)

    @property
    def low(self) -> float:
        return float(np.min(self.values))

    @property
    def high(self) -> float:
        return float(np.max(self.values))

    def sample(self) -> float:
        return float(np.random.choice(self.values))

    def clip(self, value: float) -> float:
        return float(np.clip(value, self.low, self.high))

    def __len__(self) -> int:
        return len(self.values)

    def __repr__(self) -> str:
        label = f" {self.name!r}" if self.name else ""
        return f"GeneRange({label} [{self.low}, {self.high}], n={len(self)})"


class ContinuousGeneRange(GeneRange):
    """A gene drawn from a continuous distribution rather than a discrete
    grid: sampled from N(mean, std^2) and optionally soft-bounded to
    [low, high] (default unbounded). Does not call GeneRange.__init__ (no
    `values` array to materialize) but implements the same interface, so it
    slots into ParameterSpace/clip/sample unchanged.
    """

    def __init__(
        self,
        mean: float = 0.0,
        std: float = 1.0,
        low: float = -np.inf,
        high: float = np.inf,
        name: str = "",
    ):
        self.mean = float(mean)
        self.std = float(std)
        self._low = float(low)
        self._high = float(high)
        self.name = name

    @property
    def low(self) -> float:
        return self._low

    @property
    def high(self) -> float:
        return self._high

    def sample(self) -> float:
        return self.clip(np.random.normal(self.mean, self.std))

    def clip(self, value: float) -> float:
        return float(np.clip(value, self._low, self._high))

    def __len__(self) -> int:
        return 0  # no discrete grid to count

    def __repr__(self) -> str:
        label = f" {self.name!r}" if self.name else ""
        return f"ContinuousGeneRange({label} N({self.mean}, {self.std}^2))"


class ParameterSpace:
    """Ordered collection of GeneRanges defining the genome layout."""

    def __init__(self, gene_ranges: List[GeneRange]):
        if not gene_ranges:
            raise ValueError("ParameterSpace requires at least one gene")
        self.gene_ranges = list(gene_ranges)

    @property
    def n_genes(self) -> int:
        return len(self.gene_ranges)

    def sample(self) -> np.ndarray:
        """Draw a full random genome."""
        return np.array([g.sample() for g in self.gene_ranges])

    def sample_gene(self, i: int) -> float:
        return self.gene_ranges[i].sample()

    def clip(self, genes: np.ndarray) -> np.ndarray:
        return np.array(
            [g.clip(v) for g, v in zip(self.gene_ranges, np.asarray(genes, dtype=float))]
        )

    def limits(self, i: int) -> tuple:
        gene = self.gene_ranges[i]
        return gene.low, gene.high

    def name_of(self, i: int) -> Optional[str]:
        return self.gene_ranges[i].name or None

    def __len__(self) -> int:
        return self.n_genes

    def __getitem__(self, i: int) -> GeneRange:
        return self.gene_ranges[i]

    def __repr__(self) -> str:
        return f"ParameterSpace(n_genes={self.n_genes})"
