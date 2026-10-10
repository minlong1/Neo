import os
import tempfile
import unittest

import numpy as np

from PhysicsModules.Raman.raman_neo.raman_data import RamanSpectrum, build_grid, read_spectrum
from PhysicsModules.Raman.tests.synthetic import make_spectrum, write_spectrum


class TestRamanSpectrum(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dat = write_spectrum(os.path.join(self.tmp.name, "s.dat"))

    def tearDown(self):
        self.tmp.cleanup()

    def test_build_grid(self):
        xgrid, interval = build_grid(150, 1750, 1.0)
        self.assertEqual(xgrid[0], 0.0)
        self.assertEqual(xgrid[-1], 1750.0)
        self.assertEqual(interval[0], 150)
        self.assertEqual(interval[-1], 1750)
        with self.assertRaises(ValueError):
            build_grid(500, 100, 1.0)

    def test_read_sorts_and_accepts_csv(self):
        data = make_spectrum()[::-1]  # descending x on disk
        path = write_spectrum(os.path.join(self.tmp.name, "s.csv"), data, delimiter=",")
        x, _ = read_spectrum(path)
        self.assertTrue(np.all(np.diff(x) > 0))

    def test_window_and_normalization(self):
        spec = RamanSpectrum(self.dat, xmin=150, xmax=1750, dx=1.0, normalize=True)
        self.assertEqual(spec.x[0], 150.0)
        self.assertEqual(spec.x[-1], 1750.0)
        self.assertAlmostEqual(spec.y.max(), 1.0)
        self.assertEqual(spec.get_fit_data().shape, (1601, 2))

    def test_no_normalization_keeps_intensity(self):
        spec = RamanSpectrum(self.dat, xmin=150, xmax=1750, normalize=False)
        raw = make_spectrum()
        expected = raw[raw[:, 0] == 1580.0, 1][0]
        self.assertAlmostEqual(spec.y[spec.x == 1580.0][0], expected)

    def test_baseline_subtraction(self):
        data = make_spectrum(baseline=lambda x: 0.2 + 1e-4 * x)
        path = write_spectrum(os.path.join(self.tmp.name, "bg.dat"), data)
        flat = RamanSpectrum(self.dat, xmin=150, xmax=1750, normalize=False)
        corrected = RamanSpectrum(path, xmin=150, xmax=1750, baseline_order=1, normalize=False)
        self.assertGreater(np.abs(corrected.baseline).max(), 0.1)
        self.assertLess(np.abs(corrected.y - flat.y).max(), 0.06)


if __name__ == "__main__":
    unittest.main()
