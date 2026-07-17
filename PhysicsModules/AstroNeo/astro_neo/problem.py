"""
AstroNeo as a Solvers OptimizationProblem.

This is the plug-in point between the astronomy physics (PyXspec, an X-ray
CCD spectrum, the ACX2 charge-exchange local model) and the physics-agnostic
Solvers package. Genome = model.PARAM_NAMES, the fitted model's 15 free
parameters (flat, no shared/per-path genes -- there is exactly one model).
"""

import os
import sys

from Solvers.core import GeneRange, OptimizationProblem, ParameterSpace

from PhysicsModules.AstroNeo.astro_neo.model import (
    PARAM_NAMES,
    PARAM_RANGES,
    build_xspec_model,
    genes_to_xspec_params,
)

# The original ParamsDict sampled these uniformly on the continuous interval
# [low, high] (params_type='DE' bypassed its discrete-grid branch entirely).
# Solvers.core.GeneRange represents a continuous range as a fine grid, same
# as EXAFS/NanoIndentation do for their own continuous physical parameters.
_GRID_POINTS = 100_000


def build_parameter_space() -> ParameterSpace:
    gene_ranges = []
    for name in PARAM_NAMES:
        low, high = PARAM_RANGES[name]
        step = (high - low) / _GRID_POINTS
        gene_ranges.append(GeneRange.from_bounds(low, high, step, name=name))
    return ParameterSpace(gene_ranges)


class AstroNeoProblem(OptimizationProblem):
    """Fit the model.py XSPEC spectral model to one X-ray CCD spectrum,
    scored by the Cash statistic (xspec.Fit.statistic)."""

    name = "AstroNeo"

    def __init__(
        self,
        data_dir,
        data_file,
        bg_file=None,
        rsp_file=None,
        acx2_path=None,
        xmin=7.0,
        xmax=30.0,
    ):
        self.data_dir = str(data_dir)
        self.data_file = data_file
        # None (unset) means "use whatever background/response the data
        # file's own FITS header points at" -- xspec.Spectrum()'s default,
        # passed through as its 'USE_DEFAULT' sentinel below. bg_file/
        # rsp_file only need setting to *override* that with a different
        # file than the header names.
        self.bg_file = bg_file
        self.rsp_file = rsp_file
        self.acx2_path = acx2_path
        self.xmin = xmin
        self.xmax = xmax

        self._load_spectrum()
        self.model = build_xspec_model()
        super().__init__(build_parameter_space())

    def _load_spectrum(self):
        # Imported here (not at module scope) so this module stays
        # importable without PyXspec installed, same as EXAFSProblem does
        # for larch.
        import xspec

        if self.acx2_path:
            if self.acx2_path not in sys.path:
                sys.path.append(self.acx2_path)
        import acx2_xspec  # noqa: F401 (import side effect: registers the vacx2 local model)

        xspec.xset.Xset.chatter = 0
        xspec.AllData.clear()

        old_dir = os.getcwd()
        os.chdir(self.data_dir)
        try:
            spectrum = xspec.Spectrum(
                self.data_file,
                backFile=self.bg_file or "USE_DEFAULT",
                respFile=self.rsp_file or "USE_DEFAULT",
            )
            # xLog must be set *before* .ignore(): the ignore-range string
            # below is parsed in whatever unit Plot.xAxis currently is, not
            # the unit it's later set to. Getting this order wrong silently
            # ignores every channel (0 noticed, Fit.statistic stuck at 0.0)
            # instead of raising -- caught by an end-to-end smoke run
            # against real data, not by any of the plumbing-only unit tests.
            xspec.Plot.xAxis = "angstrom"
            spectrum.ignore(f"**-{self.xmin} {self.xmax}-**")
        finally:
            os.chdir(old_dir)

        xspec.Plot.xLog = False
        xspec.Plot.yLog = False
        xspec.Plot.perHz = False
        xspec.Plot.area = True
        xspec.Plot.background = True
        xspec.Fit.statMethod = "cstat"

    def fitness(self, genes) -> float:
        # The Model is built once in __init__ and reused here -- PyXspec's
        # standard usage pattern (build once, iterate setPars(), read the
        # live-recomputed Fit.statistic). Every free parameter is set on
        # every call (genes_to_xspec_params covers all of PARAM_NAMES), so
        # there's no stale-value carryover between evaluations.
        import xspec

        self.model.setPars(genes_to_xspec_params(genes))
        return float(xspec.Fit.statistic)
