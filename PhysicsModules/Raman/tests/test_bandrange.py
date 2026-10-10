import os
import tempfile
import unittest

import numpy as np

from PhysicsModules.Raman.raman_neo.bandrange import (
    Bandrange_limits,
    build_bandrange_pars,
    check_range,
    read_bandrange_file,
    shift_grid,
)


class TestBandrange(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_defaults(self):
        br = Bandrange_limits(0)
        self.assertEqual(br.getrange_Amp()[0], 0.0)  # band can switch off
        self.assertAlmostEqual(br.getrange_Amp().max(), 1.98)
        self.assertEqual(br.getrange_Width().min(), 2.0)
        self.assertEqual(br.getrange_Cshift().min(), -15.0)
        self.assertEqual(br.get_lim()[0], 0)

    def test_custom_row(self):
        br = Bandrange_limits(1, [1, 0.1, 1.0, 0.1, 5, 10, 1, -2, 2.1, 0.5])
        np.testing.assert_allclose(br.getrange_Width(), [5, 6, 7, 8, 9])
        self.assertEqual(br.get_lim_Cshift(), (-2, 2.1, 0.5))

    def test_shift_grid(self):
        grid = shift_grid()
        self.assertEqual(len(grid), 201)
        self.assertAlmostEqual(grid[0], -10.0)
        self.assertAlmostEqual(grid[-1], 10.0)
        np.testing.assert_allclose(shift_grid((-2.0, 2.0, 1.0)), [-2, -1, 0, 1])

    def test_check_range_rejects_wide_step(self):
        with self.assertRaises(ValueError):
            check_range(0.0, 1.0, 5.0)

    def test_range_file(self):
        path = os.path.join(self.tmp.name, "ranges.csv")
        np.savetxt(
            path,
            [[0, 0.1, 1.0, 0.1, 5, 10, 1, -2, 2, 1], [1, 0.1, 2.0, 0.1, 5, 30, 1, -5, 5, 1]],
            delimiter=",",
        )
        self.assertEqual(read_bandrange_file(path, 2).shape, (2, 10))
        pars = build_bandrange_pars(2, path)
        self.assertAlmostEqual(pars[1].getrange_Width().max(), 29.0)
        with self.assertRaises(ValueError):
            read_bandrange_file(path, 3)

    def test_build_defaults(self):
        self.assertEqual(len(build_bandrange_pars(3)), 3)


if __name__ == "__main__":
    unittest.main()
