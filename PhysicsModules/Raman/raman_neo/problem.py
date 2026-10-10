"""
Raman as a Solvers OptimizationProblem.

This is the plug-in point between Raman band fitting and the physics-agnostic
Solvers package. The fit target is the (baseline-corrected, normalised)
spectrum inside the [xmin, xmax] window; the model is a sum of analytic band
line shapes (see bandshape.py). Genome layout (flat, one shared gene, like
EXAFS's shared e0):

    genes = [shift, amp_0, width_0, cshift_0, amp_1, width_1, cshift_1, ...]

``shift`` is the global wavenumber offset (the E0 analogue), shared by every
band. As in the EXAFS module it gets a dedicated grid sweep at the run's
midpoint (``on_generation_end``) and end (``on_run_end``), and once the
midpoint sweep has run, freshly sampled genomes start from the best shift
found instead of a random one (``sample_genes``).
"""

from typing import List, Optional, Sequence

import numpy as np

from Solvers.core import GeneRange, OptimizationProblem, ParameterSpace

from PhysicsModules.Raman.raman_neo.bandrange import build_bandrange_pars, shift_grid
from PhysicsModules.Raman.raman_neo.bandshape import Band, read_band_file
from PhysicsModules.Raman.raman_neo.raman_data import RamanSpectrum


def build_parameter_space(bandrange_pars: List, shift_values: Sequence[float]) -> ParameterSpace:
    """Build the flat-genome ParameterSpace from per-band Bandrange_limits."""
    gene_ranges = [GeneRange(shift_values, name="shift")]
    for i, bandrange in enumerate(bandrange_pars):
        gene_ranges.append(GeneRange(bandrange.getrange_Amp(), name=f"band{i}_amp"))
        gene_ranges.append(GeneRange(bandrange.getrange_Width(), name=f"band{i}_width"))
        gene_ranges.append(GeneRange(bandrange.getrange_Cshift(), name=f"band{i}_cshift"))
    return ParameterSpace(gene_ranges)


def band_params(genes, nband: int) -> np.ndarray:
    """Flat genome -> (nband, 4) array of ``[amp, shift, width, center_shift]``
    rows (the EXAFS ``[s02, e0, sigma2, deltaR]`` slot order, which is also
    the ``*_data.csv`` output layout)."""
    genes = np.asarray(genes, dtype=float)
    per_band = genes[1:].reshape(nband, 3)
    out = np.empty((nband, 4))
    out[:, 0] = per_band[:, 0]
    out[:, 1] = genes[0]
    out[:, 2] = per_band[:, 1]
    out[:, 3] = per_band[:, 2]
    return out


def load_bands(band_file, band_list=None, num_compounds: int = 1) -> List[Band]:
    """Instantiate the fitted Band objects from band-definition file(s).

    Single compound: ``band_file`` is one path, ``band_list`` the 1-based
    rows to fit (default: all rows), named ``Band<id>``. Multi-compound:
    ``band_file`` is a list of paths and ``band_list`` a list of row lists,
    one per compound, named ``Comp<i>Band<id>``.
    """
    bands = []
    if num_compounds > 1:
        files = _as_list(band_file)
        if len(files) != num_compounds or band_list is None or len(band_list) != num_compounds:
            raise ValueError(
                f"num_compounds={num_compounds} needs that many band files and band lists"
            )
        for i, (path, rows) in enumerate(zip(files, band_list)):
            defs = read_band_file(path)
            for band_id in rows:
                bands.append(_make_band(defs, band_id, f"Comp{i}Band{int(band_id)}"))
    else:
        defs = read_band_file(_as_list(band_file)[0])
        if band_list is None or len(band_list) == 0:
            band_list = list(range(1, len(defs) + 1))
        for band_id in band_list:
            bands.append(_make_band(defs, band_id, f"Band{int(band_id)}"))
    return bands


def _as_list(value):
    if isinstance(value, (list, tuple)):
        return [str(v).strip() for v in value]
    return [str(value).strip()]


def _make_band(defs, band_id, name) -> Band:
    idx = int(band_id) - 1
    if idx < 0 or idx >= len(defs):
        raise ValueError(f"Band {band_id} requested but the band file has {len(defs)} rows")
    d = defs[idx]
    return Band(name, d["center"], d["shape"], d["eta"], d["base_width"])


class RamanProblem(OptimizationProblem):
    """Fit a Raman spectrum as a sum of analytic bands with the NEO solvers."""

    name = "Raman"

    def __init__(
        self,
        data_file: str,
        band_file,
        band_list=None,
        xmin: float = 100.0,
        xmax: float = 1800.0,
        deltax: float = 1.0,
        weight: float = 0.0,
        baseline_order: int = 0,
        normalize: bool = True,
        bandrange_file: Optional[str] = None,
        num_compounds: int = 1,
        shift_range=(),
        shift_sweep: bool = True,
    ):
        self.data_file = data_file
        self.band_file = band_file
        self.num_compounds = int(num_compounds)
        self.weight = float(weight)
        self.shift_sweep = bool(shift_sweep)

        self.spectrum = RamanSpectrum(
            data_file,
            xmin=xmin,
            xmax=xmax,
            dx=deltax,
            baseline_order=baseline_order,
            normalize=normalize,
        )
        self.bands = load_bands(band_file, band_list, self.num_compounds)
        self.nband = len(self.bands)
        self.band_names = [band.name for band in self.bands]

        self.bandrange_pars = build_bandrange_pars(self.nband, bandrange_file)
        self.shift_values = shift_grid(shift_range)

        # Fit window, precomputed for the fitness hot path
        self.x = self.spectrum.x
        self.y = self.spectrum.y
        self.w = self.x**self.weight if self.weight != 0.0 else 1.0

        # Shift-sweep bookkeeping (EXAFS's bestE0 / secondHalf)
        self.best_shift = 0.0
        self.second_half = False

        super().__init__(build_parameter_space(self.bandrange_pars, self.shift_values))

    # ------------------------- model -------------------------

    def band_params(self, genes) -> np.ndarray:
        return band_params(genes, self.nband)

    def band_curves(self, genes, x=None) -> List[np.ndarray]:
        """Each band's curve on ``x`` (default: the fit window)."""
        x = self.x if x is None else np.asarray(x, dtype=float)
        params = self.band_params(genes)
        return [
            band.evaluate(x, amp, shift, width, cshift)
            for band, (amp, shift, width, cshift) in zip(self.bands, params)
        ]

    def model_total(self, genes, x=None) -> np.ndarray:
        """Sum of band curves on ``x`` (default: the fit window)."""
        x = self.x if x is None else np.asarray(x, dtype=float)
        genes = np.asarray(genes, dtype=float)
        shift = genes[0]
        per_band = genes[1:].reshape(self.nband, 3)
        y_total = np.zeros_like(x, dtype=float)
        for band, (amp, width, cshift) in zip(self.bands, per_band):
            if amp != 0.0:
                y_total += band.evaluate(x, amp, shift, width, cshift)
        return y_total

    def fitness(self, genes) -> float:
        """(Optionally wavenumber-weighted) sum of squared residuals over the
        fit window."""
        diff = (self.model_total(genes) - self.y) * self.w
        loss = float(np.sum(diff**2))
        if np.isnan(loss):
            return np.inf
        return loss

    def band_areas(self, genes) -> np.ndarray:
        """Integrated |intensity| of each band on the full model grid."""
        xgrid = self.spectrum.xgrid
        areas = []
        for band, row in zip(self.bands, self.band_params(genes)):
            band.set_params(*row)
            areas.append(band.area(xgrid))
        return np.asarray(areas)

    def band_optimize(self, genes, band_optimize_percent: float = 0.01) -> dict:
        """Rank bands by fractional integrated area and keep those at or above
        ``band_optimize_percent`` (the EXAFS path_optimize analogue)."""
        areas = self.band_areas(genes)
        total = areas.sum() if areas.sum() != 0 else 1.0
        fractions = areas / total
        keep = np.flatnonzero(fractions >= band_optimize_percent)
        return {
            "band_optimize_percent": band_optimize_percent,
            "original_bands": list(self.band_names),
            "optimize_bands": [self.band_names[i] for i in keep],
            "optimize_band_lists": {self.band_names[i]: float(fractions[i]) for i in keep},
        }

    # ------------------------- solver hooks -------------------------

    def sample_genes(self) -> np.ndarray:
        genes = self.space.sample()
        if self.shift_sweep and self.second_half:
            genes[0] = self.best_shift
        return genes

    def optimize_shift(self, state, population) -> float:
        """Sweep the shared shift over its grid for the global best genome,
        then apply the winner to the best and to every current individual
        (EXAFS's ``optimize_e0``). Re-sorts the population so the next
        selection sees up-to-date scores. Returns the best shift."""
        best_ind = state.globBestInd
        if best_ind is None:
            return self.best_shift

        trial = np.array(best_ind.genes, dtype=float)
        best_val = state.globBestVal
        best_shift = float(trial[0])
        for value in self.shift_values:
            trial[0] = value
            score = self.fitness(trial)
            if score < best_val:
                best_val, best_shift = score, float(value)

        self.best_shift = best_shift
        best_ind.genes[0] = best_shift
        state.globBestVal = best_val
        for individual in population.population:
            individual.genes[0] = best_shift
        population.eval_population(replace=False, sorting=True)
        return best_shift

    def on_generation_end(self, state, population) -> None:
        if self.shift_sweep and state.currGen == state.nGen // 2:
            self.optimize_shift(state, population)
            self.second_half = True

    def on_run_end(self, state, population) -> None:
        if self.shift_sweep:
            self.optimize_shift(state, population)
        # ready for another run on the same problem instance
        self.second_half = False
