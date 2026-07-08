"""
EXAFS as a Solvers OptimizationProblem.

This is the plug-in point between the EXAFS physics (larch/FEFF paths,
experimental spectrum) and the physics-agnostic Solvers package. The genome
layout is defined in individual.py: [e0, (s02, sigma2, deltaR) per path].
"""

import numpy as np

from Solvers.core import OptimizationProblem

from PhysicsModules.EXAFS.exafs_neo.individual import build_parameter_space


class _GenesView:
    """Adapts a flat gene vector to the get_e0/get_path interface that the
    EXAFS fitness function expects, without building a full Individual."""

    def __init__(self, genes):
        self.genes = np.asarray(genes, dtype=float)

    def get_e0(self):
        return float(self.genes[0])

    def get_path(self, i):
        s02, sigma2, deltaR = self.genes[1 + 3 * i : 4 + 3 * i]
        return [s02, self.get_e0(), sigma2, deltaR]


class NeoRunStateView:
    """Read-only view exposing the Solvers RunState interface over the EXAFS
    NeoPars bookkeeping, so generic operators (Metropolis, Bounded, ...) see
    live EXAFS run state without duplicated bookkeeping."""

    def __init__(self, neo_pars):
        self.neo_pars = neo_pars

    @property
    def currGen(self):
        return self.neo_pars.runPars.currGen

    @property
    def nGen(self):
        return self.neo_pars.fixedPars.nGen

    @property
    def bestDiff(self):
        return self.neo_pars.bestFitPars.bestDiff

    @property
    def currBestVal(self):
        return self.neo_pars.bestFitPars.currBestVal

    @property
    def globBestVal(self):
        return self.neo_pars.bestFitPars.globBestVal

    @property
    def diffCounter(self):
        return self.neo_pars.runPars.diffCounter

    @property
    def second_half(self):
        return self.neo_pars.runPars.secondHalf


class EXAFSProblem(OptimizationProblem):
    """The EXAFS OptimizationProblem: scores a flat genome via the larch/
    FEFF fitness in exafs_pop.py, and biases e0 sampling toward the
    best-fit value once the run passes its midpoint (`secondHalf`)."""

    name = "EXAFS"

    def __init__(self, neo_pars):
        self.neo_pars = neo_pars
        super().__init__(build_parameter_space(neo_pars.exafsRangePars.pathrange_pars))

    def fitness(self, genes) -> float:
        # Import here so the module stays importable without larch
        from PhysicsModules.EXAFS.exafs_neo.exafs_pop import fitness

        return fitness(self.neo_pars, _GenesView(genes))

    def sample_genes(self) -> np.ndarray:
        genes = self.space.sample()
        if not self.neo_pars.runPars.secondHalf:
            genes[0] = np.random.choice(self.neo_pars.exafsRangePars.rangeE0)
        else:
            genes[0] = self.neo_pars.bestFitPars.bestE0
        return genes
