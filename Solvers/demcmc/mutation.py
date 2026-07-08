"""
DE mutation operators (Forbes & Long, DE-BNN, Eq. 5-8) plus the MCMC noise
term (Eq. 48) that makes DE's selection step satisfy detailed balance when
DE is reinterpreted as an MCMC proposal/acceptance mechanism.

Notation matches the paper: for target index i, a mutant Y_i is built from
a base vector and one or two weighted difference vectors:

    rand/1:  Y_i = X_r1 + F*(X_r2 - X_r3)
    rand/2:  Y_i = X_r1 + F*(X_r2 - X_r3) + F*(X_r4 - X_r5)
    best/1:  Y_i = X_best + F*(X_r2 - X_r3)
    best/2:  Y_i = X_best + F*(X_r2 - X_r3) + F*(X_r4 - X_r5)

r1..r5 are distinct population indices, all != i.
"""

from typing import List, Sequence

import numpy as np


def _distinct_donors(n: int, exclude: int, k: int) -> np.ndarray:
    candidates = [c for c in range(n) if c != exclude]
    return np.random.choice(candidates, size=k, replace=False)


def mutate_rand_1(genomes: Sequence[np.ndarray], i: int, F: float, best: np.ndarray = None) -> np.ndarray:
    r1, r2, r3 = _distinct_donors(len(genomes), i, 3)
    return genomes[r1] + F * (genomes[r2] - genomes[r3])


def mutate_rand_2(genomes: Sequence[np.ndarray], i: int, F: float, best: np.ndarray = None) -> np.ndarray:
    r1, r2, r3, r4, r5 = _distinct_donors(len(genomes), i, 5)
    return genomes[r1] + F * (genomes[r2] - genomes[r3]) + F * (genomes[r4] - genomes[r5])


def mutate_best_1(genomes: Sequence[np.ndarray], i: int, F: float, best: np.ndarray) -> np.ndarray:
    r2, r3 = _distinct_donors(len(genomes), i, 2)
    return best + F * (genomes[r2] - genomes[r3])


def mutate_best_2(genomes: Sequence[np.ndarray], i: int, F: float, best: np.ndarray) -> np.ndarray:
    r2, r3, r4, r5 = _distinct_donors(len(genomes), i, 4)
    return best + F * (genomes[r2] - genomes[r3]) + F * (genomes[r4] - genomes[r5])


MUTATION_OPERATORS = {
    "rand/1": mutate_rand_1,
    "rand/2": mutate_rand_2,
    "best/1": mutate_best_1,
    "best/2": mutate_best_2,
}

BEST_OPERATORS = {"best/1", "best/2"}


def add_mcmc_noise(mutant: np.ndarray, sigma2: float) -> np.ndarray:
    """Eq. 48: Y_i = X_r1 + gamma*(X_r2 - X_r3) + e, e ~ N(0, sigma^2).

    Added to the mutant (regardless of which mutation operator built it) so
    the DE step satisfies the detailed-balance condition MCMC requires;
    sigma2 is kept small (paper default 1e-4) so it perturbs without
    swamping the DE differential signal.
    """
    return mutant + np.random.normal(0.0, np.sqrt(sigma2), size=mutant.shape)
