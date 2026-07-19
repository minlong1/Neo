"""
Covariance Matrix Adaptation Evolution Strategy (Hansen & Ostermeier),
the standard (mu/mu_w, lambda)-CMA-ES — not the active/mirrored-sampling
variant, which needs extra negative-weight correction terms this keeps
out of scope. Formulas follow Hansen, "The CMA Evolution Strategy: A
Tutorial" (arXiv:1604.00772); equation numbers below refer to that paper.

Each generation samples lambda=nPops candidates from a multivariate
normal N(m, sigma^2 * C), evaluates them, and uses the best mu to shift
the mean, adapt the step size sigma (cumulative step-size adaptation),
and adapt the covariance matrix C (rank-one + rank-mu update) — C is what
lets the search learn that different genes need different step scales,
rather than using one global step size for every gene.

Options (with defaults) understood on top of BaseSolver's nPops/nGen:
    sigma0 (1.0) - initial scalar step size

C is initialized as diag(((high-low)/4)**2) per gene (empirical
per-gene std of the initial population as a fallback for any gene with
infinite bounds, e.g. a default ContinuousGeneRange) rather than plain
identity. This bakes each gene's own natural scale into the search from
generation zero: a flat, gene-agnostic step size is exactly the failure
mode diagnosed in Solvers.demcmc's DE-BNN on a heterogeneous genome (its
fixed MCMC noise was ~71% of one gene's entire valid range but ~1% of
another's) — CMA-ES's whole point is to adapt this per-direction, so it
should start from a sane per-gene scale rather than learn it from an
uninformative isotropic prior.
"""

import numpy as np

from Solvers.core.base_solver import BaseSolver

MIN_POPULATION = 4  # need mu = nPops // 2 >= 2 for meaningful recombination


class CMAESSolver(BaseSolver):
    """BaseSolver wrapper around (mu/mu_w, lambda)-CMA-ES."""

    name = "CMA_ES"

    def __init__(self, problem, options=None, logger=None):
        super().__init__(problem, options, logger)
        if self.n_pops < MIN_POPULATION:
            raise ValueError(f"CMA-ES requires nPops >= {MIN_POPULATION}")

        self.sigma0 = float(self.options.get("sigma0", 1.0))

        n = self.problem.space.n_genes
        lam = self.n_pops
        mu = lam // 2
        self.n_genes = n
        self.mu = mu

        i = np.arange(1, mu + 1)
        raw_weights = np.log(mu + 0.5) - np.log(i)
        self.weights = raw_weights / raw_weights.sum()
        self.mu_eff = 1.0 / np.sum(self.weights**2)

        self.c_sigma = (self.mu_eff + 2) / (n + self.mu_eff + 5)
        self.d_sigma = (
            1
            + 2 * max(0.0, np.sqrt((self.mu_eff - 1) / (n + 1)) - 1)
            + self.c_sigma
        )
        self.c_c = (4 + self.mu_eff / n) / (n + 4 + 2 * self.mu_eff / n)
        self.c1 = 2 / ((n + 1.3) ** 2 + self.mu_eff)
        self.c_mu = min(
            1 - self.c1,
            2 * (self.mu_eff - 2 + 1 / self.mu_eff) / ((n + 2) ** 2 + self.mu_eff),
        )
        self.chi_n = np.sqrt(n) * (1 - 1 / (4 * n) + 1 / (21 * n**2))

        self.m = None
        self.sigma = None
        self.C = None
        self.p_sigma = None
        self.p_c = None

    def initialize(self) -> None:
        super().initialize()
        space = self.problem.space
        n = self.n_genes

        genes_matrix = np.array([ind.genes for ind in self.population.population])
        self.m = genes_matrix.mean(axis=0)

        scales = np.empty(n)
        for i in range(n):
            low, high = space.limits(i)
            span = high - low
            if np.isfinite(span):
                scales[i] = max(span / 4.0, 1e-8)
            else:
                std = genes_matrix[:, i].std()
                scales[i] = std if std > 0 else 1.0
        self.C = np.diag(scales**2)

        self.sigma = self.sigma0
        self.p_sigma = np.zeros(n)
        self.p_c = np.zeros(n)

    def step(self) -> None:
        space = self.problem.space
        n = self.n_genes
        lam = self.n_pops

        C_sym = (self.C + self.C.T) / 2
        eigvals, B = np.linalg.eigh(C_sym)
        eigvals = np.clip(eigvals, 1e-20, None)
        D = np.sqrt(eigvals)

        Z = np.random.normal(size=(lam, n))
        Y = Z @ np.diag(D) @ B.T
        X = self.m + self.sigma * Y
        X_clipped = np.array([space.clip(x) for x in X])

        for i, individual in enumerate(self.population.population):
            individual.set(X_clipped[i])
        scores = self.population.eval_population()

        y_used = (X_clipped - self.m) / self.sigma

        order = np.argsort(scores)
        best_idx = order[: self.mu]
        y_best = y_used[best_idx]
        y_w = self.weights @ y_best

        self.m = self.m + self.sigma * y_w

        c_inv_sqrt = B @ np.diag(1.0 / D) @ B.T
        self.p_sigma = (1 - self.c_sigma) * self.p_sigma + np.sqrt(
            self.c_sigma * (2 - self.c_sigma) * self.mu_eff
        ) * (c_inv_sqrt @ y_w)

        self.sigma *= np.exp(
            (self.c_sigma / self.d_sigma) * (np.linalg.norm(self.p_sigma) / self.chi_n - 1)
        )

        gen = self.state.currGen
        h_sigma_lhs = np.linalg.norm(self.p_sigma) / np.sqrt(
            1 - (1 - self.c_sigma) ** (2 * gen)
        )
        h_sigma = 1.0 if h_sigma_lhs < (1.4 + 2 / (n + 1)) * self.chi_n else 0.0

        self.p_c = (1 - self.c_c) * self.p_c + h_sigma * np.sqrt(
            self.c_c * (2 - self.c_c) * self.mu_eff
        ) * y_w

        rank_mu_update = sum(
            self.weights[j] * np.outer(y_best[j], y_best[j]) for j in range(self.mu)
        )
        delta_h_sigma = (1 - h_sigma) * self.c_c * (2 - self.c_c)
        self.C = (
            (1 - self.c1 - self.c_mu) * self.C
            + self.c1 * (np.outer(self.p_c, self.p_c) + delta_h_sigma * self.C)
            + self.c_mu * rank_mu_update
        )
        self.C = (self.C + self.C.T) / 2
