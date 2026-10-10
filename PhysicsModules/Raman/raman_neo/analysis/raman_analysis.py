"""
Post-fit analysis and reporting.

Raman analogue of the EXAFS module's ``analysis/EXAFS_Analysis.py``: turns a
finished fit (a RamanProblem plus its best genome) into a peak table, a
fit/residual curve, summary statistics, and an optional plot.

    from PhysicsModules.Raman.raman_neo.analysis.raman_analysis import RamanAnalysis
    analysis = RamanAnalysis(problem, result.best_individual.genes)
    print(analysis.report())
"""

import csv

import numpy as np

from PhysicsModules.Raman.raman_neo.analysis.raman_score import r_factor

PEAK_TABLE_FIELDS = [
    "band",
    "shape",
    "nominal_center",
    "fitted_center",
    "amplitude",
    "width",
    "center_shift",
    "shift",
    "area",
    "area_fraction",
]


class RamanAnalysis:
    """Summarise one genome of a RamanProblem (usually the best fit)."""

    def __init__(self, problem, genes, fitness=None):
        self.problem = problem
        self.genes = np.asarray(genes, dtype=float)
        self.params = problem.band_params(self.genes)  # (nband, 4)
        self.shift = float(self.genes[0])
        self.fitness = problem.fitness(self.genes) if fitness is None else float(fitness)

        spectrum = problem.spectrum
        self.xgrid = spectrum.xgrid
        self.exp = spectrum.exp
        self.interval = spectrum.interval
        self.model = problem.model_total(self.genes, x=self.xgrid)

    @classmethod
    def from_result(cls, problem, result):
        """Build from a Solvers SolverResult."""
        return cls(problem, result.best_individual.genes, fitness=result.best_value)

    def peak_table(self):
        """Per-band dicts with fitted peak descriptors."""
        areas = self.problem.band_areas(self.genes)
        total = areas.sum() if areas.sum() != 0 else 1.0
        rows = []
        for band, (amp, shift, width, cshift), area in zip(self.problem.bands, self.params, areas):
            rows.append(
                {
                    "band": band.name,
                    "shape": band.shape,
                    "nominal_center": band.center,
                    "fitted_center": band.center + cshift + shift,
                    "amplitude": float(amp),
                    "width": float(width),
                    "center_shift": float(cshift),
                    "shift": float(shift),
                    "area": float(area),
                    "area_fraction": float(area / total),
                }
            )
        return rows

    def write_peak_table(self, path) -> None:
        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=PEAK_TABLE_FIELDS)
            writer.writeheader()
            writer.writerows(self.peak_table())

    def r_factor(self):
        return r_factor(self.model, self.exp, self.interval)

    def residual(self):
        return self.exp - self.model

    def report(self):
        """Formatted multi-line text report."""
        lines = [
            "=" * 66,
            "Raman Neo fit report",
            "=" * 66,
            f"Global shift (cm^-1) : {self.shift:.3f}",
            f"Final fitness        : {self.fitness:.6g}",
            f"R-factor             : {self.r_factor():.5f}",
            "-" * 66,
            f"{'band':>12} {'center':>9} {'amp':>8} {'width':>8} {'area':>10} {'frac':>7}",
        ]
        for row in self.peak_table():
            lines.append(
                f"{row['band']:>12} {row['fitted_center']:>9.2f} "
                f"{row['amplitude']:>8.3f} {row['width']:>8.2f} "
                f"{row['area']:>10.2f} {100 * row['area_fraction']:>6.1f}%"
            )
        lines.append("=" * 66)
        return "\n".join(lines)

    def plot(self, ax=None, show_bands=False):
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots(figsize=(8, 5))
        idx = self.interval
        x = self.xgrid[idx]
        ax.plot(x, self.exp[idx], "k-", lw=1.2, label="experiment")
        ax.plot(x, self.model[idx], "r--", lw=1.2, label="fit")
        if show_bands:
            for band, curve in zip(self.problem.bands, self.problem.band_curves(self.genes, x=x)):
                ax.plot(x, curve, lw=0.8, alpha=0.7, label=band.name)
        offset = 0.2 * (np.max(self.exp[idx]) if idx.size else 1.0)
        ax.plot(x, self.residual()[idx] - offset, "b-", lw=0.8, label="residual")
        ax.set_xlabel("Raman shift (cm$^{-1}$)")
        ax.set_ylabel("Intensity (norm.)" if self.problem.spectrum.normalize else "Intensity")
        ax.legend()
        return ax
