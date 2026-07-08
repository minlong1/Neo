"""
Differential Evolution solver (DE/rand/1/bin), Storn & Price 1997.

The algorithm itself lives in `differential_evolution_step`, decoupled from
BaseSolver's own generation loop so a physics module that runs its own loop
(e.g. EXAFS's ExafsNeo.run()) can call it per-generation directly against
its own population, rather than duplicating the DE math. See
PhysicsModules/EXAFS/exafs_neo/neoSolver.py:NeoSolver_DE for that wiring.

Options (with defaults) understood on top of BaseSolver's nPops/nGen:
    F (0.5)   - differential (mutation) weight, typically in [0, 2]
    CR (0.9)  - crossover probability, in [0, 1]

For each target vector x_i, a trial is built from three other, distinct
population members a, b, c as:

    mutant = a + F * (b - c)
    trial_j = mutant_j if rand() < CR or j == j_rand else x_i_j

(j_rand guarantees at least one gene comes from the mutant, so trial is
never a plain copy of the target). The trial is clipped back into the
parameter space's bounds, then replaces the target if it scores better
(one-to-one greedy selection, evaluated against a snapshot of the
population so replacements within a generation don't affect each other's
donor vectors).
"""

import numpy as np

from Solvers.core.base_solver import BaseSolver

MIN_POPULATION = 4  # 1 target + 3 distinct donors


def differential_evolution_step(population, F: float = 0.5, CR: float = 0.9) -> None:
    """Advance any Population-like object by one DE/rand/1/bin generation.

    `population` needs only: `.problem` (`.space.clip`, `.fitness`),
    `.population` (list of individuals exposing `.genes` and `len()`),
    `.eval_population(replace, sorting)` returning per-individual scores in
    population order, and `.generate_individual()` to construct a new,
    correctly-typed individual (its genes are immediately overwritten via
    `.set(...)`, so what generate_individual() samples doesn't matter).
    Both Solvers.core.Population and EXAFS's NeoPopulations satisfy this.
    """
    problem = population.problem
    space = problem.space

    target_scores = population.eval_population(replace=False, sorting=False)
    targets = list(population.population)
    n = len(targets)
    if n < MIN_POPULATION:
        raise ValueError(
            f"DE requires at least {MIN_POPULATION} individuals "
            f"(3 distinct donors per target), got {n}"
        )

    next_population = []
    for i, target in enumerate(targets):
        donors = np.random.choice(
            [c for c in range(n) if c != i], size=3, replace=False
        )
        a, b, c = (targets[d] for d in donors)
        mutant = a.genes + F * (b.genes - c.genes)

        n_genes = len(target)
        trial_genes = target.genes.copy()
        j_rand = np.random.randint(n_genes)
        cross_mask = np.random.random(n_genes) < CR
        cross_mask[j_rand] = True
        trial_genes = np.where(cross_mask, mutant, trial_genes)
        trial_genes = space.clip(trial_genes)

        trial_score = problem.fitness(trial_genes)
        if trial_score < target_scores[i]:
            trial = population.generate_individual()
            trial.set(trial_genes)
            next_population.append(trial)
        else:
            next_population.append(target)

    population.population = next_population
    population.eval_population()


class DESolver(BaseSolver):
    """BaseSolver wrapper around `differential_evolution_step` — for
    physics modules that hand control to a solver's own run() loop rather
    than driving DE per-generation themselves."""

    name = "DE"

    def __init__(self, problem, options=None, logger=None):
        super().__init__(problem, options, logger)
        if self.n_pops < MIN_POPULATION:
            raise ValueError(
                f"DE requires nPops >= {MIN_POPULATION} (3 distinct donors per target)"
            )
        self.F = float(self.options.get("F", 0.5))
        self.CR = float(self.options.get("CR", 0.9))

    def step(self):
        differential_evolution_step(self.population, self.F, self.CR)
