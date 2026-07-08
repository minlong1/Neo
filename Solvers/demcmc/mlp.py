"""
Minimal numpy-only feedforward MLP: the "neural network" side of DE-NN /
DE-BNN (Forbes & Long, "DE-BNN: An evolutionary approach to Bayesian
neural network posterior sampling", Neurocomputing 678 (2026) 133103).

Matches the paper's matrix formulation (their Eq. 26-32) for an arbitrary
number of layers: for layer sizes (n_0, n_1, ..., n_L),

    S^k = Z^{k-1} @ W_{k-1} + b_{k-1}      (Z^0 = X)
    Z^k = activation(S^k)   for k < L, output_activation(S^k) for k == L

Weight matrices and bias vectors are kept as separate, individually
shaped arrays (weights[k]: n_k x n_{k+1}, biases[k]: n_{k+1}) — the "matrix
structure" the paper's SVD refinement needs (Section 3.5) — but are
flattened to/from a single gene vector for DE's mutation/crossover, which
operate elementwise on flat arrays.
"""

from typing import List, Sequence, Tuple

import numpy as np

from Solvers.core.parameter_space import ContinuousGeneRange, ParameterSpace


def relu(x: np.ndarray) -> np.ndarray:
    return np.maximum(0.0, x)


def identity(x: np.ndarray) -> np.ndarray:
    return x


class MLPStructure:
    """Layer-size layout plus flatten/unflatten/forward for a plain MLP.

    layer_sizes: (n_in, n_hidden_1, ..., n_hidden_k, n_out).
    """

    def __init__(
        self,
        layer_sizes: Sequence[int],
        activation=relu,
        output_activation=identity,
    ):
        self.layer_sizes = list(layer_sizes)
        if len(self.layer_sizes) < 2:
            raise ValueError("layer_sizes needs at least (n_in, n_out)")
        self.activation = activation
        self.output_activation = output_activation

        # Each entry: (kind, shape, slice-into-the-flat-gene-vector)
        self._layout: List[Tuple[str, Tuple[int, ...], slice]] = []
        offset = 0
        for n_in, n_out in zip(self.layer_sizes[:-1], self.layer_sizes[1:]):
            w_size = n_in * n_out
            self._layout.append(("W", (n_in, n_out), slice(offset, offset + w_size)))
            offset += w_size
            self._layout.append(("b", (n_out,), slice(offset, offset + n_out)))
            offset += n_out
        self.n_genes = offset

    @property
    def n_layers(self) -> int:
        return len(self.layer_sizes) - 1

    def unflatten(self, genes: np.ndarray) -> Tuple[List[np.ndarray], List[np.ndarray]]:
        """genes (flat, length n_genes) -> (weights, biases), matrix-shaped."""
        genes = np.asarray(genes, dtype=float)
        weights, biases = [], []
        for kind, shape, sl in self._layout:
            arr = genes[sl].reshape(shape)
            (weights if kind == "W" else biases).append(arr)
        return weights, biases

    def flatten(self, weights: Sequence[np.ndarray], biases: Sequence[np.ndarray]) -> np.ndarray:
        """Inverse of unflatten: matrix-shaped (weights, biases) -> flat genes."""
        parts = []
        wi = bi = 0
        for kind, _shape, _sl in self._layout:
            if kind == "W":
                parts.append(np.asarray(weights[wi], dtype=float).ravel())
                wi += 1
            else:
                parts.append(np.asarray(biases[bi], dtype=float).ravel())
                bi += 1
        return np.concatenate(parts)

    def forward(self, genes: np.ndarray, X: np.ndarray) -> np.ndarray:
        """X: (n_samples, n_in) -> y_hat: (n_samples, n_out)."""
        weights, biases = self.unflatten(genes)
        Z = np.asarray(X, dtype=float)
        last = len(weights) - 1
        for k, (W, b) in enumerate(zip(weights, biases)):
            S = Z @ W + b
            Z = self.output_activation(S) if k == last else self.activation(S)
        return Z

    def build_parameter_space(self, weight_std: float = None, bias_std: float = 0.1) -> ParameterSpace:
        """He-normal initialization (paper Section 6.2: "weight matrices and
        bias vector are a type of He initialization"): each weight fan-in
        gets std = sqrt(2 / n_in); biases default to a small fixed std.
        `weight_std`, if given, overrides He-init with a fixed std for
        every weight gene instead.
        """
        gene_ranges = []
        for kind, shape, _sl in self._layout:
            if kind == "W":
                n_in = shape[0]
                std = weight_std if weight_std is not None else np.sqrt(2.0 / n_in)
                gene_ranges.extend(
                    ContinuousGeneRange(std=std) for _ in range(shape[0] * shape[1])
                )
            else:
                gene_ranges.extend(ContinuousGeneRange(std=bias_std) for _ in range(shape[0]))
        return ParameterSpace(gene_ranges)
