import csv
import os
import tempfile
import unittest

import numpy as np

from PhysicsModules.Raman.raman_neo.analysis import raman_score as rs
from PhysicsModules.Raman.raman_neo.analysis.raman_analysis import RamanAnalysis
from PhysicsModules.Raman.raman_neo.problem import RamanProblem
from PhysicsModules.Raman.tests.synthetic import (
    TRUE_BANDS,
    true_genes,
    write_band_file,
    write_spectrum,
)

HERE = os.path.dirname(os.path.abspath(__file__))
SAMPLE_DATA = os.path.join(os.path.dirname(HERE), "raman_neo", "data")


class AnalysisTestBase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.dat = write_spectrum(os.path.join(cls.tmp.name, "spectrum.dat"))
        cls.bands = write_band_file(os.path.join(cls.tmp.name, "bands.txt"))
        cls.problem = RamanProblem(
            cls.dat, cls.bands, band_list=[1, 2, 3, 4, 5, 6], xmin=150, xmax=1750, normalize=False
        )

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()


class TestRamanAnalysis(AnalysisTestBase):
    def test_truth_report(self):
        analysis = RamanAnalysis(self.problem, true_genes(6))
        self.assertLess(analysis.fitness, 1e-12)
        self.assertLess(analysis.r_factor(), 1e-12)
        np.testing.assert_allclose(analysis.residual()[analysis.interval], 0.0, atol=1e-9)

        table = analysis.peak_table()
        self.assertEqual(len(table), 6)
        for row, (center, amp, width, _, _) in zip(table, TRUE_BANDS):
            self.assertEqual(row["fitted_center"], center)
            self.assertAlmostEqual(row["amplitude"], amp)
            self.assertAlmostEqual(row["width"], width)
        self.assertAlmostEqual(sum(r["area_fraction"] for r in table), 1.0)
        self.assertEqual(table[4]["area"], 0.0)  # decoy

        report = analysis.report()
        self.assertIn("Raman Neo fit report", report)
        self.assertIn("Band4", report)

    def test_write_peak_table(self):
        path = os.path.join(self.tmp.name, "peaks.csv")
        RamanAnalysis(self.problem, true_genes(6)).write_peak_table(path)
        with open(path) as f:
            rows = list(csv.DictReader(f))
        self.assertEqual(len(rows), 6)
        self.assertEqual(float(rows[3]["fitted_center"]), 1580.0)

    def test_from_result(self):
        class FakeResult:
            class best_individual:
                genes = true_genes(6)

            best_value = 0.25

        analysis = RamanAnalysis.from_result(self.problem, FakeResult)
        self.assertEqual(analysis.fitness, 0.25)

    def test_plot(self):
        try:
            import matplotlib

            matplotlib.use("Agg")
        except ImportError:
            self.skipTest("matplotlib not installed")
        ax = RamanAnalysis(self.problem, true_genes(6)).plot(show_bands=True)
        self.assertGreaterEqual(len(ax.lines), 3 + 6)


class TestRamanScore(AnalysisTestBase):
    def test_build_grid(self):
        xg, iv = rs.build_grid(150, 1750, 1.0)
        self.assertEqual(xg[0], 0.0)
        self.assertEqual(xg[-1], 1750.0)
        self.assertEqual(iv[0], 150)
        self.assertEqual(iv[-1], 1750)

    def test_standalone_score_matches_problem_fitness(self):
        xg, iv = rs.build_grid(150, 1750, 1.0)
        exp = rs.load_experimental(self.dat, xg, normalize=False, interval=iv)
        bands = rs.bands_from_file(self.bands, [1, 2, 3, 4, 5, 6])
        genes = true_genes(6)
        genes[1] = 0.3  # make it imperfect
        model = rs.synthesize(bands, self.problem.band_params(genes), xg)
        self.assertAlmostEqual(rs.score(model, exp, xg, iv), self.problem.fitness(genes))

    def test_bundled_sample_truth_fits(self):
        xg, iv = rs.build_grid(150, 1750, 1.0)
        exp = rs.load_experimental(
            os.path.join(SAMPLE_DATA, "sample_raman.dat"), xg, baseline_order=1, normalize=True, interval=iv
        )
        bands = rs.bands_from_file(os.path.join(SAMPLE_DATA, "sample_bands.txt"), [1, 2, 3, 4])
        truth = np.array([[0.45, 0, 9, 0], [0.30, 0, 12, 0], [0.80, 0, 28, 0], [1.00, 0, 18, 0]])
        model = rs.synthesize(bands, truth, xg)
        self.assertLess(rs.r_factor(model, exp, iv), 0.1)

    def test_score_zero_for_identical(self):
        y = np.array([0.0, 1.0, 2.0, 3.0])
        self.assertEqual(rs.score(y, y), 0.0)
        self.assertEqual(rs.r_factor(y, y), 0.0)


if __name__ == "__main__":
    unittest.main()
