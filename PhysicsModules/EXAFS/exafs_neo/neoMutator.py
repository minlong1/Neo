"""
EXAFS-facing mutator wrapper. The mutation algorithms live in
Solvers.ga.mutators; this module adapts them to the historical EXAFS API.

The EXAFS-specific e0 mutation during the second half of the run stays
here — it acts on the shared-e0 genome convention, which is EXAFS physics,
not solver logic.
"""

import numpy as np

from Solvers.ga.mutators import create_mutator

# Historical (mutOpt attribute, mutType label) per mut_options value.
# Note the mutOpt attribute has always been mut_options + 1 (see the
# original per-class mutOpt assignments in EXAFS Neo <= 2.x).
_MUT_LABELS = {
    0: (1, "Mutate Per Individual"),
    1: (2, "Mutate Per Path"),
    2: (3, "Mutate Per Trait"),
    3: (4, "Mutate Metropolis"),
    4: (5, "Mutate Bounded Per Range"),
    5: (6, "Mutate DE"),
}


class _DEMutatorStub:
    """Placeholder for DE mutation (was ExafsMutator_DE, never implemented).

    DE now lives in Solvers.de as a full solver rather than a GA mutator.
    """

    def __init__(self, mut_chance, logger=None):
        self.mut_chance = mut_chance
        self.logger = logger

    def mutate(self, pops):
        pass

    def __str__(self):
        return f"mutation chance: {self.mut_chance}"


class NeoMutator:
    def __init__(self, logger=None):
        self.mutator = None
        self.logger = logger
        self.mutator_type = None
        self.exafs_pars = None
        self.mutator_score = 0  # TODO: need to check if this is needed

    def initialize(self, exafs_pars):
        self.exafs_pars = exafs_pars

        self.mutator_type = exafs_pars.mutPars.mutOpt
        if self.mutator_type not in _MUT_LABELS:
            raise ValueError("Invalid mutator type")

        if self.mutator_type == 5:
            self.mutator = _DEMutatorStub(
                exafs_pars.mutPars.mutChance, logger=self.logger
            )
        else:
            self.mutator = create_mutator(
                self.mutator_type,
                mut_chance=exafs_pars.mutPars.mutChance,
                logger=self.logger,
            )
        # Preserve the historical attribute names on the operator
        self.mutator.mutOpt, self.mutator.mutType = _MUT_LABELS[self.mutator_type]
        self.mutator.mutChance = exafs_pars.mutPars.mutChance
        self.mutator.mutChanceE0 = exafs_pars.mutPars.mutChanceE0

        return self.mutator

    def __str__(self):
        if self.mutator is None:
            return "None Mutator selected"
        else:
            return f"Mutator Type: {self.mutator.mutType}, {self.mutator}"

    def mutate(self, pops):
        if self.mutator is None:
            raise ValueError("Mutator is not initialized")
        else:
            if self.exafs_pars.runPars.secondHalf:
                self.mutate_e0(pops)
            self.mutator.mutate(pops)

    def mutate_e0(self, pops):
        """
        Mutate the e0 value in the second half
        """
        if np.random.random() * 100 < self.mutator.mutChanceE0:
            e0 = np.random.choice(self.exafs_pars.exafsRangePars.rangeE0)
            if self.exafs_pars.verbose_lvl >= 5:
                self.logger.print(f"Mutate e0 to: {e0:.3f}")
            for individual in pops.population:
                individual.set_e0(e0)
