"""
Posterior sample collection and prediction (Forbes & Long, DE-BNN,
Section 5). DE-MCMC treats the DE population as (a set of) Markov chains:
at each generation past burn-in, the accepted candidate(s) at one or more
indices are recorded as draws from the posterior over weights/biases.

Single-chain (num_chains=1) records only the population-best candidate
each generation (Section 5.2); multi-chain records the top `num_chains`
candidates, pooling all of them as posterior samples — DE's population-
based nature means these are already meaningfully different chains
because each is a combination of differentials from the others.
"""

from typing import List

import numpy as np


class PosteriorResult:
    """Accumulated posterior samples (flat genomes) collected past
    burn-in, plus prediction/summary methods over them."""

    def __init__(self):
        self.samples: List[np.ndarray] = []

    def collect(self, population, num_chains: int = 1) -> None:
        top = population.population_sorted[:num_chains]
        for individual, _score in top:
            self.samples.append(individual.genes.copy())

    def as_array(self) -> np.ndarray:
        if not self.samples:
            raise ValueError("No posterior samples collected yet (past burn_in?)")
        return np.array(self.samples)

    def mean_genes(self) -> np.ndarray:
        """Elementwise posterior mean of each weight/bias component."""
        return self.as_array().mean(axis=0)

    def mode_genes(self, bins: int = 50) -> np.ndarray:
        """Elementwise posterior mode, approximated per-component via a
        histogram (the paper's "mode prediction", Section 6.2.1)."""
        arr = self.as_array()
        modes = np.empty(arr.shape[1])
        for j in range(arr.shape[1]):
            hist, edges = np.histogram(arr[:, j], bins=bins)
            k = int(np.argmax(hist))
            modes[j] = 0.5 * (edges[k] + edges[k + 1])
        return modes

    def predict(self, problem, X=None):
        """Predictive posterior mean (Eq. 49-50): forward-pass every
        posterior sample and average the resulting predictions — not the
        prediction of the averaged weights. Returns (mean_prediction,
        per_sample_predictions) so callers can also derive credible
        intervals from the same forward passes.
        """
        preds = np.array([problem.predict(genes, X) for genes in self.samples])
        return preds.mean(axis=0), preds

    def credible_interval(self, problem, X=None, low: float = 5.0, high: float = 95.0):
        """Percentile interval of the predictive posterior at each input
        (e.g. low=5, high=95 for a 5%-95% band, matching the paper's
        figures)."""
        _, preds = self.predict(problem, X)
        return np.percentile(preds, low, axis=0), np.percentile(preds, high, axis=0)

    def __len__(self) -> int:
        return len(self.samples)
