"""
NanoIndentation as a Solvers OptimizationProblem.

This is the plug-in point between nanoindentation mechanics and the
physics-agnostic Solvers package. The fit target is the sliced unloading
segment of a load-displacement curve; the model is a sum of Oliver-Pharr
power laws y = A(h - hf)^m. Genome layout (flat, no shared genes):

    genes = [A_0, hf_0, m_0, A_1, hf_1, m_1, ...]   (3 genes per path)
"""

import numpy as np

from Solvers.core import GeneRange, OptimizationProblem, ParameterSpace

from PhysicsModules.NanoIndentation.nanoindentation_neo.nano_neo_data import (
    NanoIndent_Data,
)
from PhysicsModules.NanoIndentation.nanoindentation_neo.pathObj import (
    DEFAULT_A_RANGE,
    DEFAULT_HF_RANGE,
    DEFAULT_M_RANGE,
    get_model,
    select_opt,
)


def build_parameter_space(npaths, A_range=(), hf_range=(), m_range=()) -> ParameterSpace:
    """Build the flat-genome ParameterSpace: (A, hf, m) grids per path.

    Ranges are (low, high, step) triples; empty ranges fall back to the
    historical nano_indent defaults.
    """
    A_values = select_opt(A_range, DEFAULT_A_RANGE)
    hf_values = select_opt(hf_range, DEFAULT_HF_RANGE)
    m_values = select_opt(m_range, DEFAULT_M_RANGE)

    gene_ranges = []
    for i in range(npaths):
        gene_ranges.append(GeneRange(A_values, name=f"path{i}_A"))
        gene_ranges.append(GeneRange(hf_values, name=f"path{i}_hf"))
        gene_ranges.append(GeneRange(m_values, name=f"path{i}_m"))
    return ParameterSpace(gene_ranges)


class NanoIndentationProblem(OptimizationProblem):
    """Fit nanoindentation unloading curves with the NEO solvers."""

    name = "NanoIndentation"

    def __init__(
        self,
        data_file: str,
        data_cutoff=(0.1, 0.9),
        npaths: int = 1,
        A_range=(),
        hf_range=(),
        m_range=(),
        fits: str = "OliverPharr",
    ):
        self.npaths = npaths
        self.model = get_model(fits)

        self.data_file = data_file
        self.data_obj = NanoIndent_Data(data_file)
        self.data_obj.pre_processing(limits=(data_cutoff[0], data_cutoff[1]))
        self.x_slice = self.data_obj.get_slice_data()[:, 0]
        self.y_slice = self.data_obj.get_slice_data()[:, 1]

        super().__init__(build_parameter_space(npaths, A_range, hf_range, m_range))

    def model_total(self, genes) -> np.ndarray:
        """Sum of per-path model curves over the fitted x slice."""
        genes = np.asarray(genes, dtype=float).reshape(self.npaths, 3)
        y_total = np.zeros_like(self.x_slice)
        with np.errstate(invalid="ignore"):
            for A, h_f, m in genes:
                y_total = y_total + self.model.func(self.x_slice, A, h_f, m)
        return y_total

    def fitness(self, genes) -> float:
        """Sum of squared residuals against the sliced unloading data.

        (h - hf)^m is NaN where h < hf with fractional m; such genomes are
        scored inf so sorting discards them (the original implementation
        dropped NaN individuals and regenerated instead).
        """
        y_total = self.model_total(genes)
        loss = np.sum((y_total - self.y_slice) ** 2)
        if np.isnan(loss):
            return np.inf
        return float(loss)
