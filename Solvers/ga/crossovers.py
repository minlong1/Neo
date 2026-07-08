"""
GA crossover operators, ported from EXAFS neoCrossOver.py, operating on
flat gene vectors instead of per-path tuples.

Numeric option IDs preserved: 0 = uniform, 1 = single point, 2 = dual point,
3 = arithmetic, 4 = or, 5 = average.
"""

import numpy as np

from Solvers.core.individual import Individual
from Solvers.core.population import Population


class CrossoverBase:
    """Breeds pairs drawn from `next_population` to refill the population
    back to size; subclasses implement `crossover_pair()`, the per-gene
    combination rule."""

    crossover_type = None
    croType = "Base"

    def __init__(self, n_cross: int, logger=None):
        self.logger = logger
        self.n_cross = n_cross

    def crossover_pair(self, pops: Population, individual1: Individual, individual2: Individual) -> Individual:
        raise NotImplementedError

    def crossover(self, pops: Population) -> None:
        """Fill the population back up by breeding pairs from next_population."""
        if len(pops.next_population) <= 2:
            return
        children = []
        for _ in range(self.n_cross):
            par_ind = np.random.choice(len(pops.next_population), size=2, replace=False)
            ind1 = pops.next_population[par_ind[0]]
            ind2 = pops.next_population[par_ind[1]]
            children.append(self.crossover_pair(pops, ind1, ind2))

        pops.next_population.extend(children)
        pops.population = pops.next_population

    def __str__(self):
        return f"Crossover Option: {self.croType}"


class UniformCrossover(CrossoverBase):
    crossover_type = 0
    croType = "Uniform Crossover"

    def crossover_pair(self, pops, individual1, individual2):
        mask = np.random.randint(0, 2, size=len(individual1)).astype(bool)
        genes = np.where(mask, individual1.genes, individual2.genes)
        return Individual(individual1.space, genes)


class SinglePointCrossover(CrossoverBase):
    crossover_type = 1
    croType = "Single Point Crossover"

    def crossover_pair(self, pops, individual1, individual2, co_point: int = None):
        n = len(individual1)
        if co_point is None:
            co_point = np.random.randint(1, n) if n > 1 else 0
        co_point = int(np.clip(co_point, 0, n))
        genes = np.concatenate(
            [individual1.genes[:co_point], individual2.genes[co_point:]]
        )
        return Individual(individual1.space, genes)


class DualPointCrossover(CrossoverBase):
    # TODO: implement dual point crossover (stub carried over from neoCrossOver.py)
    crossover_type = 2
    croType = "Dual Point Crossover"

    def crossover_pair(self, pops, individual1, individual2):
        pass


class ArithmeticCrossover(CrossoverBase):
    crossover_type = 3
    croType = "Arithmetic Crossover"

    def crossover_pair(self, pops, individual1, individual2):
        mask = np.logical_and(
            np.random.randint(0, 2, size=len(individual1)),
            np.random.randint(0, 2, size=len(individual1)),
        )
        genes = np.where(mask, individual1.genes, individual2.genes)
        return Individual(individual1.space, genes)


class OrCrossover(CrossoverBase):
    crossover_type = 4
    croType = "Or Crossover"

    def crossover_pair(self, pops, individual1, individual2):
        mask = np.logical_or(
            np.random.randint(0, 2, size=len(individual1)),
            np.random.randint(0, 2, size=len(individual1)),
        )
        genes = np.where(mask, individual1.genes, individual2.genes)
        return Individual(individual1.space, genes)


class AverageCrossover(CrossoverBase):
    crossover_type = 5
    croType = "Average Crossover"

    def crossover_pair(self, pops, individual1, individual2):
        # TODO check if averaged values go out of bounds (carried over)
        genes = (individual1.genes + individual2.genes) / 2
        return Individual(individual1.space, genes)


_CROSSOVERS = {
    0: UniformCrossover,
    1: SinglePointCrossover,
    2: DualPointCrossover,
    3: ArithmeticCrossover,
    4: OrCrossover,
    5: AverageCrossover,
}


def create_crossover(cro_opt: int, n_cross: int, logger=None) -> CrossoverBase:
    if cro_opt not in _CROSSOVERS:
        raise ValueError(f"Invalid crossover type: {cro_opt}")
    return _CROSSOVERS[cro_opt](n_cross, logger=logger)
