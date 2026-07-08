"""
Bayesian neural network regression as a Solvers OptimizationProblem — the
plug-in point between DE-BNN's DEMCMCSolver and a training dataset.

The genome is the flattened weight/bias vector of an MLPStructure; fitness
is sum-of-squared-error, matching the paper's DE fitness (a Gaussian
likelihood assumption is what makes SSE minimization equivalent to MAP/
posterior-mode estimation in the Bayesian framing of Section 5).
"""

import numpy as np

from Solvers.core.problem import OptimizationProblem
from Solvers.demcmc.mlp import MLPStructure, identity, relu


class BNNRegressionProblem(OptimizationProblem):
    """Reference DE-BNN problem: fits an MLPStructure's flattened weights
    to (X, y) by minimizing SSE. The reference OptimizationProblem for
    DEMCMCSolver."""

    name = "BNN"

    def __init__(
        self,
        layer_sizes,
        X,
        y,
        activation=relu,
        output_activation=identity,
        weight_std: float = None,
        bias_std: float = 0.1,
    ):
        self.mlp = MLPStructure(layer_sizes, activation=activation, output_activation=output_activation)
        self.X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float)
        n_out = layer_sizes[-1]
        self.y = y.reshape(-1, n_out) if y.ndim == 1 and n_out == 1 else y.reshape(-1, n_out)

        super().__init__(self.mlp.build_parameter_space(weight_std=weight_std, bias_std=bias_std))

    def fitness(self, genes: np.ndarray) -> float:
        y_pred = self.mlp.forward(genes, self.X)
        return float(np.sum((y_pred - self.y) ** 2))

    def predict(self, genes: np.ndarray, X=None) -> np.ndarray:
        return self.mlp.forward(genes, self.X if X is None else np.asarray(X, dtype=float))
