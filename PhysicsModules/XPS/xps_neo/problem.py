"""
XPS as a Solvers OptimizationProblem — scaffold.

This module is the plug-in point between XPS physics (X-ray Photoelectron
Spectroscopy peak fitting) and the physics-agnostic Solvers package. It is
not implemented yet; see the module README for the contract to fill in.
"""

import numpy as np

from Solvers.core import GeneRange, OptimizationProblem, ParameterSpace


class XPSProblem(OptimizationProblem):
    """Scaffold for fitting XPS spectra with the NEO solvers.

    To implement:
    1. Build the ParameterSpace: one GeneRange per fitting parameter
       (e.g. per peak: binding energy, amplitude, FWHM, Gaussian/Lorentzian
       mix), each with the discrete grid of allowed values.
    2. Load the experimental spectrum in __init__ (keep all I/O here — the
       solver never sees files).
    3. Implement fitness(genes): synthesize the model spectrum from the gene
       vector and return a scalar misfit against the measured spectrum
       (lower is better).
    4. Optionally override sample_genes() to bias initialization, and the
       on_generation_end/on_run_end hooks to write per-generation outputs.

    See PhysicsModules/EXAFS/exafs_neo/problem.py for a complete example.
    """

    name = "XPS"

    def __init__(self):
        # Placeholder single-gene space so the scaffold is instantiable;
        # replace with the real fitting-parameter grid.
        super().__init__(ParameterSpace([GeneRange([0.0], name="placeholder")]))

    def fitness(self, genes: np.ndarray) -> float:
        raise NotImplementedError("XPS fitting is not implemented yet")
