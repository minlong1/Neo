"""
EXAFS-facing crossover wrapper. The crossover algorithms live in
Solvers.ga.crossovers and operate on flat gene vectors; this module adapts
them to the historical EXAFS API and rebuilds full EXAFS individuals
(with path metadata) from the generic children.
"""

import numpy as np

from Solvers.ga.crossovers import create_crossover

# Historical croType labels, indexed by croOpt
_CRO_TYPES = {
    0: "Uniform Crossover",
    1: "Single Point Crossover",
    2: "Dual Point Crossover",
    3: "Arithmetic Crossover",
    4: "Or Crossover",
    5: "Average Crossover",
}


class NeoCrossover:
    """Back-compat shim preserving the historical `NeoCrossover().initialize(exafs_pars)`
    / `.crossover(pops)` API; delegates the actual gene-mixing to a
    Solvers.ga crossover operator selected by `croOpt`."""

    def __init__(self, logger=None):
        self.logger = logger
        self.exafs_pars = None
        self.crossover_type = None
        self.crossover_operator = None
        self.crossover_score = 0  # TODO: maybe implement this?

    def initialize(self, exafs_pars):
        self.exafs_pars = exafs_pars

        self.crossover_type = exafs_pars.crossPars.croOpt
        self.crossover_operator = create_crossover(
            self.crossover_type,
            n_cross=exafs_pars.selPars.nCross,
            logger=self.logger,
        )
        # Preserve the historical attribute names on the operator
        self.crossover_operator.croOpt = self.crossover_type
        self.crossover_operator.croType = _CRO_TYPES[self.crossover_type]

    def __str__(self):
        if self.crossover_operator is None:
            return "Crossover is not selected"
        else:
            return f"Crossover Type: {self.crossover_type}, {self.crossover_operator}"

    def crossover(self, pops):
        if self.crossover_operator is None:
            raise ValueError("Crossover is not initialized")
        else:
            temp_population = []
            if len(pops.next_population) > 2:
                for _ in range(self.exafs_pars.selPars.nCross):
                    par_ind = np.random.choice(
                        len(pops.next_population), size=2, replace=False
                    )
                    ind1 = pops.next_population[par_ind[0]]
                    ind2 = pops.next_population[par_ind[1]]
                    child = self.crossover_single(pops, ind1, ind2)
                    temp_population.append(child)

                pops.next_population.extend(temp_population)
                pops.population = pops.next_population

    def crossover_single(self, pops, ind1, ind2):
        if self.crossover_operator is None:
            raise ValueError("Crossover is not initialized")
        else:
            proto = self.crossover_operator.crossover_pair(pops, ind1, ind2)
            child = pops.generate_individual()
            child.set(proto.genes)
            return child
