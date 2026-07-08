"""
GA mutation operators, ported from EXAFS neoMutator.py.

Numeric option IDs preserved from NeoMutator.initialize:
0 = per individual, 1 = per gene (was "per path"), 2 = per trait (stub),
3 = Metropolis, 4 = bounded.

All chances are probabilities in [0, 1].
"""

import copy

import numpy as np

from Solvers.core.population import Population


class MutatorBase:
    """Perturbs a population in place after crossover; subclasses implement
    `mutate()`. `mut_chance` is a probability in [0, 1] (not the
    percent-scaled convention some EXAFS mutation paths historically use)."""

    mutator_type = None
    mutType = "Base"

    def __init__(self, mut_chance: float = 0.3, logger=None):
        self.logger = logger
        self.mut_chance = mut_chance
        self.nmut = 0

    def mutate(self, pops: Population):
        pass

    def __str__(self):
        return f"mutation chance: {self.mut_chance}"


class PerIndividualMutator(MutatorBase):
    """Replace whole individuals with fresh random ones."""

    mutator_type = 0
    mutType = "Mutate Per Individual"

    def mutate(self, pops: Population):
        for i, _ in enumerate(pops.population):
            if np.random.random() < self.mut_chance:
                pops.population[i] = pops.generate_individual()
                self.nmut += 1


class PerGeneMutator(MutatorBase):
    """Resample individual genes (generalizes EXAFS 'mutate per path')."""

    mutator_type = 1
    mutType = "Mutate Per Gene"

    def mutate(self, pops: Population):
        for i, individual in enumerate(pops.population):
            if np.random.random() < self.mut_chance:
                individual.mutate(self.mut_chance)
                self.nmut += 1


class PerTraitMutator(MutatorBase):
    # TODO: implement per-trait mutation (stub carried over from neoMutator.py)
    mutator_type = 2
    mutType = "Mutate Per Trait"

    def mutate(self, pops: Population):
        pass


class MetropolisMutator(MutatorBase):
    """Accept mutations that improve fitness, or worsen it with a
    temperature-controlled probability (simulated-annealing style)."""

    mutator_type = 3
    mutType = "Mutate Metropolis"

    def mutate(self, pops: Population):
        state = pops.state
        for i, indi in enumerate(pops.population):
            if np.random.random() < self.mut_chance:
                og_indi = copy.deepcopy(indi)
                og_score = pops.problem.fitness(og_indi.genes)
                mut_indi = copy.deepcopy(indi)
                mut_indi.mutate(self.mut_chance)
                mut_score = pops.problem.fitness(mut_indi.genes)

                T = -state.bestDiff / np.log(1 - (state.currGen / state.nGen))
                if mut_score < og_score:
                    new_indi = mut_indi
                elif np.exp(-(mut_score - og_score) / T) > np.random.uniform():
                    new_indi = mut_indi
                else:
                    new_indi = og_indi

                pops.population[i] = new_indi


class BoundedMutator(MutatorBase):
    """Shift genes toward their bounds with a generation-shrinking step
    (non-uniform mutation)."""

    mutator_type = 4
    mutType = "Mutate Bounded Per Range"

    def mutate(self, pops: Population):
        state = pops.state
        space = pops.problem.space

        def delta_fun(t, delta_val):
            rnd = np.random.random()
            return delta_val * (1 - rnd ** (1 - (t / state.nGen)) ** 5)

        for indi in pops.population:
            if np.random.random() < self.mut_chance:
                for j in range(len(indi)):
                    low, high = space.limits(j)
                    val = indi.get_gene(j)
                    if np.random.randint(2) == 0:
                        val = val + delta_fun(state.currGen, high - val)
                    else:
                        val = val - delta_fun(state.currGen, val - low)
                    indi.set_gene(j, val)


_MUTATORS = {
    0: PerIndividualMutator,
    1: PerGeneMutator,
    2: PerTraitMutator,
    3: MetropolisMutator,
    4: BoundedMutator,
}


def create_mutator(mut_opt: int, mut_chance: float = 0.3, logger=None) -> MutatorBase:
    if mut_opt not in _MUTATORS:
        raise ValueError(f"Invalid mutator type: {mut_opt}")
    return _MUTATORS[mut_opt](mut_chance, logger=logger)
