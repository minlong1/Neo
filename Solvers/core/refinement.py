"""
Genome-agnostic periodic local-search refinement: perturb a candidate's
genes with small per-gene noise, keep the perturbation only if it beats
the *population's* current best fitness (not just the candidate's own).

Generalizes Solvers.demcmc.refinement.local_search_refine, which does the
same thing but round-trips through an MLPStructure's unflatten/flatten
purely to reshape the flat genome into per-layer weight matrices before
perturbing — a step that's a no-op for the perturbation itself (adding
noise to a reshaped-then-flattened vector is identical to adding it to
the flat vector directly), so this version works against any
OptimizationProblem, not just BNN-shaped ones. Solvers.demcmc.refinement
is untouched; this is a new, separate, additive function.

`perturb_scale` is a *fraction of each gene's own (high - low) range*,
not an absolute noise width — computed once per gene from
`problem.space.limits(i)` (falling back to a fixed absolute width for any
gene with infinite bounds, e.g. a default ContinuousGeneRange). A flat,
gene-agnostic noise width is exactly the failure mode diagnosed in
Solvers.demcmc's DE-BNN against a heterogeneous genome (its fixed MCMC
noise scale was ~71% of one gene's entire valid range but only ~1% of
another's) — this makes the function scale-aware by construction instead
of repeating that mistake.
"""

import numpy as np


def _gene_perturb_widths(space, perturb_scale: float, fallback_width: float = 1.0) -> np.ndarray:
    n = space.n_genes
    widths = np.empty(n)
    for i in range(n):
        low, high = space.limits(i)
        span = high - low
        widths[i] = span * perturb_scale if np.isfinite(span) else fallback_width
    return widths


def local_search_refine(
    population,
    problem,
    population_min: float,
    individuals=None,
    n_samples: int = 20,
    perturb_scale: float = 0.05,
) -> int:
    """Perturb each of `individuals` (default: the whole population)
    `n_samples` times with uniform noise scaled to each gene's own
    range, keeping the best perturbation only if it beats
    `population_min`. Returns the number of individuals updated.
    """
    space = problem.space
    targets = population.population if individuals is None else individuals
    widths = _gene_perturb_widths(space, perturb_scale)
    n_updated = 0

    for individual in targets:
        best_genes = individual.genes
        best_score = problem.fitness(individual.genes)
        improved = False

        for _ in range(n_samples):
            noise = np.random.uniform(-widths, widths)
            trial_genes = space.clip(individual.genes + noise)
            score = problem.fitness(trial_genes)
            if score < population_min and score < best_score:
                best_score = score
                best_genes = trial_genes
                improved = True

        if improved:
            individual.set(best_genes)
            n_updated += 1

    if n_updated:
        population.eval_population()
    return n_updated
