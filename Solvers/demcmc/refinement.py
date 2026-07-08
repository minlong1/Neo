"""
Periodic refinement techniques (Forbes & Long, DE-BNN, Section 3.5-3.7),
run every `refine_every` generations rather than every generation (the
paper's own reasoning: "reduce computational expense... at the same time
preserve the general trend of the search", Section 3.8).

Both implemented techniques accept a modified candidate only if it beats
the population's current best fitness (not just the individual's own
fitness) — this is the paper's stated criterion for both SVD ("lower than
the population minimum", Section 3.5) and local search ("less than the
population best fitness", Section 3.7), a deliberately strict bar meant to
inject rare, truly-better candidates rather than accept lateral moves.

Clustering (Section 3.6: k-means/spectral/agglomerative refinement) is not
implemented — it would pull in a clustering dependency (e.g. scikit-learn)
this package doesn't otherwise need, for one of three interchangeable
refinement techniques. `cluster_refine` below is a documented stub.
"""

import numpy as np

# Eq. 18-20 discrete parameter grids for the three SVD variants.
LAMBDA_CHOICES = np.round(np.arange(0.0, 2.0, 0.1), 2)  # scalar multiply
KAPPA_CHOICES = np.round(np.arange(1.05, 1.091, 0.01), 3)  # exponent
ETA_CHOICES = np.round(np.arange(4.0, 5.01, 0.05), 3)  # log-then-exponent


def _svd_variants(W: np.ndarray) -> list:
    """The three modified reconstructions of W from Eq. 18-20, each
    obtained by transforming W's singular values then rebuilding W."""
    if min(W.shape) == 0:
        return []
    U, S, Vt = np.linalg.svd(W, full_matrices=False)

    lam = np.random.choice(LAMBDA_CHOICES)
    W1 = U @ np.diag(lam * S) @ Vt  # Eq. 18

    kappa = np.random.choice(KAPPA_CHOICES)
    W2 = U @ np.diag(S**kappa) @ Vt  # Eq. 19

    eta = np.random.choice(ETA_CHOICES)
    W3 = U @ np.diag(np.log(S + 1.0) ** eta) @ Vt  # Eq. 20

    return [W1, W2, W3]


def svd_refine(population, problem, population_min: float) -> int:
    """Try each SVD variant on each weight matrix of each individual;
    keep the modification if it beats `population_min`. Returns the
    number of individuals updated.

    Requires `problem.mlp` (an MLPStructure) for flatten/unflatten —
    i.e. this is specific to BNNRegressionProblem-shaped problems.
    """
    space = problem.space
    mlp = problem.mlp
    n_updated = 0

    for individual in population.population:
        weights, biases = mlp.unflatten(individual.genes)
        best_genes = individual.genes
        best_score = problem.fitness(individual.genes)
        improved = False

        for wi, W in enumerate(weights):
            for W_mod in _svd_variants(W):
                trial_weights = list(weights)
                trial_weights[wi] = W_mod
                trial_genes = space.clip(mlp.flatten(trial_weights, biases))
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


def local_search_refine(
    population,
    problem,
    population_min: float,
    n_samples: int = 20,
    iteration: int = 0,
    perturb_scale: float = 0.1,
) -> int:
    """Perturb each individual's weight matrices with uniform noise
    n_samples times (Eq. 24), keeping the best perturbation if it beats
    `population_min`. n_samples grows with the length of the current
    stagnation streak (Eq. 25: n_samples * (iteration/1000 + 1)).
    Returns the number of individuals updated.
    """
    space = problem.space
    mlp = problem.mlp
    n_samples_now = int(n_samples * (iteration / 1000.0 + 1))
    n_updated = 0

    for individual in population.population:
        weights, biases = mlp.unflatten(individual.genes)
        best_genes = individual.genes
        best_score = problem.fitness(individual.genes)
        improved = False

        for _ in range(n_samples_now):
            trial_weights = [
                W + np.random.uniform(-perturb_scale, perturb_scale, size=W.shape)
                for W in weights
            ]
            trial_genes = space.clip(mlp.flatten(trial_weights, biases))
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


def cluster_refine(population, problem, population_min: float, clustering_type: str = "kmeans"):
    """Clustering-based refinement (Section 3.6): cluster each population
    of weight-matrix/bias-vector candidates (k-means, spectral, or
    agglomerative), combine per-cluster centers into proposed candidates,
    and keep any that improve on the population minimum.

    Not implemented — see module docstring. Raises NotImplementedError so
    callers fail loudly rather than silently no-op.
    """
    raise NotImplementedError(
        "cluster_refine (DE-BNN Section 3.6) is not implemented; use "
        "svd_refine and/or local_search_refine, both of which are."
    )
