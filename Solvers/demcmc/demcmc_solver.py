"""
DE-MCMC / DE-BNN solver (Forbes & Long, "DE-BNN: An evolutionary approach
to Bayesian neural network posterior sampling", Neurocomputing 678 (2026)
133103): differential evolution reinterpreted as a Markov Chain Monte
Carlo sampler, used here to train a neural network's weights/biases while
simultaneously producing posterior samples for Bayesian prediction.

DE's own selection step (accept the trial if it scores better) already
matches an MCMC acceptance rule; with small Gaussian noise added to the
mutant for detailed balance (mutation.add_mcmc_noise, Eq. 48), the
sequence of accepted candidates at each population index forms a Markov
chain, and the population as a whole gives NP such chains "for free".

Options (with defaults) on top of BaseSolver's nPops/nGen:
    F (0.5), CR (0.7), mutation_operator ("rand/1" | "rand/2" | "best/1" |
        "best/2"), sigma2 (1e-4)              - base DE-MCMC parameters
    burn_in (0)                               - generations discarded
                                                before posterior collection
    num_chains (1)                            - top-N candidates per
                                                generation kept as samples
    hyper_mutation (False), stagnation_window (10)
                                               - resample F/CR/operator
                                                (Section 3.1-3.4) while the
                                                running-average residual
                                                is ~0
    svd (False), local_search (False), refine_every (10),
    local_search_samples (20)                 - periodic refinement
                                                (Section 3.5, 3.7)
"""

import numpy as np

from Solvers.core.base_solver import BaseSolver
from Solvers.demcmc.hyper_mutation import StagnationTracker, sample_hyperparameters
from Solvers.demcmc.mutation import BEST_OPERATORS, MUTATION_OPERATORS, add_mcmc_noise
from Solvers.demcmc.posterior import PosteriorResult
from Solvers.demcmc.refinement import local_search_refine, svd_refine

# rand/2 and best/2 need 5 distinct donors besides the target; since
# hyper-mutation can switch operators mid-run, require enough population
# for the most demanding operator regardless of the configured default.
MIN_POPULATION = 6


def de_mcmc_generation_step(population, F: float, CR: float, mutation_operator: str, sigma2: float = 1e-4) -> None:
    """One DE-MCMC generation against any Population-like object (see
    Solvers.de.differential_evolution_step for the exact interface
    required: .problem, .population, .eval_population(replace, sorting),
    .generate_individual()).
    """
    if mutation_operator not in MUTATION_OPERATORS:
        raise ValueError(
            f"Unknown mutation_operator {mutation_operator!r}; "
            f"available: {sorted(MUTATION_OPERATORS)}"
        )
    operator_fn = MUTATION_OPERATORS[mutation_operator]
    needs_best = mutation_operator in BEST_OPERATORS

    problem = population.problem
    space = problem.space

    target_scores = population.eval_population(replace=False, sorting=False)
    targets = list(population.population)
    genomes = [t.genes for t in targets]
    n = len(targets)
    if n < MIN_POPULATION:
        raise ValueError(
            f"DE-MCMC requires at least {MIN_POPULATION} individuals "
            f"(rand/2 and best/2 need 5 distinct donors), got {n}"
        )

    best_genome = genomes[int(np.argmin(target_scores))] if needs_best else None

    next_population = []
    for i, target in enumerate(targets):
        mutant = operator_fn(genomes, i, F, best_genome)
        mutant = add_mcmc_noise(mutant, sigma2)

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


class DEMCMCSolver(BaseSolver):
    """Ties one DE-MCMC generation, periodic SVD/local-search refinement,
    and post-burn-in posterior collection into a single BaseSolver.
    Requires nPops >= MIN_POPULATION (6) since hyper-mutation can switch
    to best/2 or rand/2 mid-run regardless of the configured default
    operator, and both need 5 distinct donors."""

    name = "DE_MCMC"

    def __init__(self, problem, options=None, logger=None):
        super().__init__(problem, options, logger)
        if self.n_pops < MIN_POPULATION:
            raise ValueError(f"DE-MCMC requires nPops >= {MIN_POPULATION}")

        self.F = float(self.options.get("F", 0.5))
        self.CR = float(self.options.get("CR", 0.7))
        self.mutation_operator = self.options.get("mutation_operator", "rand/1")
        if self.mutation_operator not in MUTATION_OPERATORS:
            raise ValueError(f"Unknown mutation_operator {self.mutation_operator!r}")
        self.sigma2 = float(self.options.get("sigma2", 1e-4))

        self.burn_in = int(self.options.get("burn_in", 0))
        self.num_chains = int(self.options.get("num_chains", 1))

        self.hyper_mutation = bool(self.options.get("hyper_mutation", False))
        self._stagnation = StagnationTracker(window=int(self.options.get("stagnation_window", 10)))

        self.svd = bool(self.options.get("svd", False))
        self.local_search = bool(self.options.get("local_search", False))
        self.refine_every = int(self.options.get("refine_every", 10))
        self.local_search_samples = int(self.options.get("local_search_samples", 20))

        self.posterior = PosteriorResult()

    def step(self):
        F, CR, operator = self.F, self.CR, self.mutation_operator
        if self.hyper_mutation and self._stagnation.is_stagnant():
            F, CR, operator = sample_hyperparameters()

        de_mcmc_generation_step(self.population, F, CR, operator, sigma2=self.sigma2)
        self._stagnation.update(self.population.population_sorted[0][1])

        if (self.svd or self.local_search) and self.state.currGen % self.refine_every == 0:
            population_min = self.population.population_sorted[0][1]
            if self.svd:
                svd_refine(self.population, self.problem, population_min)
            if self.local_search:
                local_search_refine(
                    self.population,
                    self.problem,
                    self.population.population_sorted[0][1],
                    n_samples=self.local_search_samples,
                    iteration=self._stagnation.iteration,
                )

        if self.state.currGen > self.burn_in:
            self.posterior.collect(self.population, self.num_chains)
