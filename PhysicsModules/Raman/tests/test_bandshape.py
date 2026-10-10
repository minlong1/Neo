import os
import tempfile
import unittest

import numpy as np

from PhysicsModules.Raman.raman_neo.bandshape import (
    Band,
    evaluate_shape,
    gaussian,
    lorentzian,
    polynomial_baseline,
    pseudo_voigt,
    read_band_file,
)


class TestLineShapes(unittest.TestCase):
    def setUp(self):
        self.x = np.linspace(0.0, 200.0, 2001)

    def test_unit_peak_height(self):
        for func in (lorentzian, gaussian, pseudo_voigt):
            with self.subTest(func=func.__name__):
                y = func(self.x, 100.0, 5.0)
                self.assertAlmostEqual(y.max(), 1.0, places=6)
                self.assertAlmostEqual(self.x[np.argmax(y)], 100.0)

    def test_lorentzian_hwhm(self):
        self.assertAlmostEqual(lorentzian(np.array([105.0]), 100.0, 5.0)[0], 0.5)

    def test_gaussian_sigma(self):
        self.assertAlmostEqual(gaussian(np.array([105.0]), 100.0, 5.0)[0], np.exp(-0.5))

    def test_voigt_eta_limits(self):
        np.testing.assert_allclose(pseudo_voigt(self.x, 100, 5, 1.0), lorentzian(self.x, 100, 5))
        np.testing.assert_allclose(pseudo_voigt(self.x, 100, 5, 0.0), gaussian(self.x, 100, 5))

    def test_unknown_shape_rejected(self):
        with self.assertRaises(ValueError):
            evaluate_shape("square", self.x, 100.0, 5.0)
        with self.assertRaises(ValueError):
            Band("b", 100.0, shape="square")


class TestBand(unittest.TestCase):
    def test_absolute_center_and_area(self):
        x = np.arange(0.0, 2000.0, 0.5)
        band = Band("b", 500.0, "lorentzian")
        band.set_params(amp=2.0, shift=3.0, width=4.0, center_shift=-1.0)
        self.assertEqual(band.absolute_center(), 502.0)
        y = band.band2chi(x)
        self.assertAlmostEqual(x[np.argmax(y)], 502.0)
        # Lorentzian area = amp * pi * HWHM (minus truncated tails)
        self.assertAlmostEqual(band.area(x), 2.0 * np.pi * 4.0, delta=0.1)

    def test_evaluate_is_stateless(self):
        band = Band("b", 500.0)
        x = np.array([500.0])
        self.assertEqual(band.evaluate(x, 1.5, 0.0, 5.0, 0.0)[0], 1.5)
        self.assertEqual(band.amp, 0.0)


class TestBandFile(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, text):
        path = os.path.join(self.tmp.name, "bands.txt")
        with open(path, "w") as fh:
            fh.write(text)
        return path

    def test_defaults_and_comments(self):
        path = self.write("# header\n\n520.0\n1350, voigt, 0.7, 25\n")
        bands = read_band_file(path)
        self.assertEqual(len(bands), 2)
        self.assertEqual(
            bands[0], {"center": 520.0, "shape": "lorentzian", "eta": 0.5, "base_width": 10.0}
        )
        self.assertEqual(bands[1]["shape"], "voigt")
        self.assertEqual(bands[1]["eta"], 0.7)

    def test_empty_file_rejected(self):
        with self.assertRaises(ValueError):
            read_band_file(self.write("# nothing\n"))


class TestBaseline(unittest.TestCase):
    def test_zero_order_is_zero(self):
        x = np.linspace(0, 10, 11)
        np.testing.assert_array_equal(polynomial_baseline(x, x, order=0), np.zeros_like(x))

    def test_recovers_linear_background_under_peak(self):
        x = np.linspace(100.0, 1800.0, 1701)
        background = 0.1 + 1e-4 * x
        y = background + lorentzian(x, 900.0, 10.0)
        baseline = polynomial_baseline(x, y, order=1)
        self.assertLess(np.max(np.abs(baseline - background)), 0.02)


if __name__ == "__main__":
    unittest.main()
