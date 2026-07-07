import os
import tempfile
import unittest

import numpy as np

from PhysicsModules.NanoIndentation.nanoindentation_neo.nano_neo_data import (
    NanoIndent_Data,
)
from PhysicsModules.NanoIndentation.tests.synthetic import make_curve, write_csv


class TestNanoIndentData(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.csv = write_csv(os.path.join(self.tmp.name, "curve.csv"))

    def tearDown(self):
        self.tmp.cleanup()

    def test_csv_fallback_type(self):
        data = NanoIndent_Data(self.csv)
        self.assertEqual(data._ftype, "NAN")

    def test_raw_data_shape(self):
        data = NanoIndent_Data(self.csv)
        data.read_data()
        self.assertEqual(data.get_raw_data().shape, make_curve().shape)

    def test_preprocessing_slices_unloading_segment(self):
        data = NanoIndent_Data(self.csv)
        data.pre_processing(limits=(0.1, 0.9))
        sliced = data.get_slice_data()

        raw = make_curve()
        peak_x = raw[np.argmax(raw[:, 0]), 0]
        unload_y = raw[np.argmax(raw[:, 0]) :, 1]

        # slice must lie in the unloading branch, strictly inside the limits
        self.assertGreater(len(sliced), 0)
        self.assertLessEqual(sliced[:, 0].max(), peak_x)
        delta = abs(unload_y[0] - unload_y[-1])
        self.assertGreater(sliced[:, 1].min(), unload_y[-1] + 0.1 * delta - 1e-9)
        self.assertLess(sliced[:, 1].max(), unload_y[0] - 0.1 * delta + 1e-9)

    def test_tighter_limits_give_smaller_slice(self):
        wide = NanoIndent_Data(self.csv)
        wide.pre_processing(limits=(0.05, 0.95))
        tight = NanoIndent_Data(self.csv)
        tight.pre_processing(limits=(0.3, 0.7))
        self.assertLess(tight.get_length(), wide.get_length())


if __name__ == "__main__":
    unittest.main()
