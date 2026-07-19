"""
Particle Swarm Optimization (Kennedy & Eberhart 1995), global-best variant
with Clerc & Kennedy's (2002) constriction coefficients.

Options (with defaults) understood on top of BaseSolver's nPops/nGen:
    w (0.7298)   - inertia weight
    c1 (1.49618) - cognitive coefficient (pull toward the individual's own
                   best-known position)
    c2 (1.49618) - social coefficient (pull toward the swarm's best-known
                   position)
    local_search (False), refine_every (10), local_search_samples (10)
                 - optional periodic local-search refinement of the
                   current best individual (Solvers.core.refinement.
                   local_search_refine), targeting PSO's actual weakness:
                   not its mean/best, but run-to-run variance (it can land
                   in a mediocre basin and stay there). Refines only
                   population.best() every `refine_every` generations, so
                   the extra cost stays small (n_samples extra evals per
                   refine_every generations, not per whole population).

Each individual i tracks a velocity v_i alongside its position x_i (=
Individual.genes), its own best-known position p_i, and the swarm shares a
single global-best position g. Every generation:

    v_i = w*v_i + c1*r1*(p_i - x_i) + c2*r2*(g - x_i)   # r1, r2 ~ U(0,1), per-gene
    x_i = clip(x_i + v_i)

then re-evaluate and update p_i/g. The (w, c1, c2) defaults are the
standard constriction values, chosen so the swarm converges without a
separate velocity-clamp (vMax) parameter.
"""

import numpy as np

from Solvers.core.base_solver import BaseSolver
from Solvers.core.refinement import local_search_refine

MIN_POPULATION = 2  # a swarm of 1 has no distinct personal/global best


class PSOSolver(BaseSolver):
    """BaseSolver wrapper around global-best PSO."""

    name = "PSO"

    def __init__(self, problem, options=None, logger=None):
        super().__init__(problem, options, logger)
        if self.n_pops < MIN_POPULATION:
            raise ValueError(f"PSO requires nPops >= {MIN_POPULATION}")
        self.w = float(self.options.get("w", 0.7298))
        self.c1 = float(self.options.get("c1", 1.49618))
        self.c2 = float(self.options.get("c2", 1.49618))

        self.local_search = bool(self.options.get("local_search", False))
        self.refine_every = int(self.options.get("refine_every", 10))
        self.local_search_samples = int(self.options.get("local_search_samples", 10))

        self.velocities = None
        self.personal_best_genes = None
        self.personal_best_scores = None

    def initialize(self) -> None:
        super().initialize()
        n_genes = self.problem.space.n_genes
        self.velocities = np.zeros((self.n_pops, n_genes))

        scores = self.population.eval_population(replace=False, sorting=False)
        self.personal_best_genes = np.array(
            [individual.genes.copy() for individual in self.population.population]
        )
        self.personal_best_scores = np.array(scores, dtype=float)

    def step(self) -> None:
        space = self.problem.space
        scores = self.population.eval_population(replace=False, sorting=False)

        for i, score in enumerate(scores):
            if score < self.personal_best_scores[i]:
                self.personal_best_scores[i] = score
                self.personal_best_genes[i] = self.population.population[i].genes.copy()

        g_best_idx = int(np.argmin(self.personal_best_scores))
        g_best_genes = self.personal_best_genes[g_best_idx]

        positions = np.array([individual.genes for individual in self.population.population])
        n_genes = space.n_genes
        r1 = np.random.random((self.n_pops, n_genes))
        r2 = np.random.random((self.n_pops, n_genes))

        self.velocities = (
            self.w * self.velocities
            + self.c1 * r1 * (self.personal_best_genes - positions)
            + self.c2 * r2 * (g_best_genes - positions)
        )
        new_positions = positions + self.velocities

        for i, individual in enumerate(self.population.population):
            individual.set(space.clip(new_positions[i]))

        self.population.eval_population()

        if self.local_search and self.state.currGen % self.refine_every == 0:
            population_min = self.population.population_sorted[0][1]
            local_search_refine(
                self.population,
                self.problem,
                population_min,
                individuals=[self.population.best()],
                n_samples=self.local_search_samples,
            )
