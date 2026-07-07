"""
Generic individual: a flat gene vector tied to a ParameterSpace.
"""

import numpy as np

from Solvers.core.parameter_space import ParameterSpace


class Individual:
    def __init__(self, space: ParameterSpace, genes: np.ndarray = None):
        self.space = space
        if genes is None:
            self.genes = space.sample()
        else:
            genes = np.asarray(genes, dtype=float)
            if genes.shape != (space.n_genes,):
                raise ValueError(
                    f"Expected {space.n_genes} genes, got shape {genes.shape}"
                )
            self.genes = genes.copy()

    def get(self) -> np.ndarray:
        return self.genes.copy()

    def set(self, genes: np.ndarray) -> None:
        genes = np.asarray(genes, dtype=float)
        if genes.shape != self.genes.shape:
            raise ValueError(f"Expected shape {self.genes.shape}, got {genes.shape}")
        self.genes = genes.copy()

    def get_gene(self, i: int) -> float:
        return float(self.genes[i])

    def set_gene(self, i: int, value: float) -> None:
        self.genes[i] = value

    def mutate(self, chance: float) -> int:
        """Resample each gene with probability `chance` (in [0, 1]).

        Returns the number of genes mutated.
        """
        n_mutated = 0
        for i in range(len(self.genes)):
            if np.random.random() < chance:
                self.genes[i] = self.space.sample_gene(i)
                n_mutated += 1
        return n_mutated

    def copy(self) -> "Individual":
        return Individual(self.space, self.genes)

    def __len__(self) -> int:
        return len(self.genes)

    def __getitem__(self, i: int) -> float:
        return float(self.genes[i])

    def __str__(self) -> str:
        return f"Individual with genes {self.genes}"
