"""
Genetic-algorithm solvers, ported from EXAFS neoSolver.py.

Options (with defaults) understood on top of BaseSolver's nPops/nGen:
    selOpt (0), croOpt (0), mutOpt (1), nBestSample (0.3), nLuckSample (0.2),
    mutChance (0.3)
"""

import numpy as np

from Solvers.core.base_solver import BaseSolver
from Solvers.ga.selectors import create_selector
from Solvers.ga.crossovers import create_crossover
from Solvers.ga.mutators import create_mutator


def rechenberg_update(curr_gen, diff_counter, glob_best_val, curr_best_val, mut_chance):
    """Rechenberg's 1/5 success rule: adapt the mutation chance based on how
    often recent generations improved the best fit.

    Returns the updated (diff_counter, mut_chance). Pure function so physics
    modules with their own bookkeeping can reuse the exact same schedule.
    """
    best_diff = np.abs(glob_best_val - curr_best_val)
    if curr_gen > 20:
        if best_diff < 0.1:
            diff_counter += 1
        else:
            diff_counter -= 1
        diff_counter = int(np.clip(diff_counter, 0, None))

        ratio = abs(diff_counter) / float(curr_gen)
        if ratio > 0.2:
            if (mut_chance + 0.0025) < 1.0:
                mut_chance = abs(mut_chance + 0.0025)
        elif ratio < 0.2:
            if (mut_chance - 0.0025) > 0:
                mut_chance -= 0.0025
        mut_chance = float(np.clip(mut_chance, 0, 1.0))

    return diff_counter, mut_chance


class GASolver(BaseSolver):
    name = "GA"

    def __init__(self, problem, options=None, logger=None):
        super().__init__(problem, options, logger)

        self.selector = create_selector(
            self.options.get("selOpt", 0),
            n_pops=self.n_pops,
            n_best_sample=self.options.get("nBestSample", 0.3),
            n_lucky_sample=self.options.get("nLuckSample", 0.2),
            logger=logger,
        )
        self.crossover = create_crossover(
            self.options.get("croOpt", 0),
            n_cross=self.selector.nCross,
            logger=logger,
        )
        self.mutator = create_mutator(
            self.options.get("mutOpt", 1),
            mut_chance=self.options.get("mutChance", 0.3),
            logger=logger,
        )

    def step(self):
        self.selector.select(self.population)
        self.crossover.crossover(self.population)
        self.mutator.mutate(self.population)
        self.population.eval_population()


class GARechenbergSolver(GASolver):
    """GA with Rechenberg's 1/5 success rule adapting the mutation chance."""

    name = "GA_Rechenberg"

    def step(self):
        self.selector.select(self.population)
        self.crossover.crossover(self.population)
        self.rechenberg_mutation()
        self.mutator.mutate(self.population)
        self.population.eval_population()

    def rechenberg_mutation(self):
        state = self.state
        state.diffCounter, self.mutator.mut_chance = rechenberg_update(
            state.currGen,
            state.diffCounter,
            state.globBestVal,
            state.currBestVal,
            self.mutator.mut_chance,
        )
