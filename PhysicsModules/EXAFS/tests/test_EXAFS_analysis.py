"""
Unit test for EXAFS Analysis

The fixture read by these tests, tests/cu_test_files/cu_results/run1_data.csv,
is real GA output (13-path fit against tests/cu_test_files/cu_paths/cu_10k.xmu,
solOpt=0, nPops=60, 14 generations completed) — it was missing from the repo
entirely (read_result_files silently failed with an UnboundLocalError on the
missing directory; see EXAFS_Analysis.py's FileNotFoundError guard). The
expected values below are pinned to that checked-in fixture (a
characterization test): regenerating the fixture from a different GA run
requires updating them to match.
"""
import sys
import larch
import os
import unittest
import numpy as np
from pathlib import Path

sys.path.append('gui')
from PhysicsModules.EXAFS.exafs_neo.analysis import EXAFS_Analysis
from PhysicsModules.EXAFS.exafs_neo.analysis import larch_score


class TestCase(unittest.TestCase):
    mylarch = larch.Interpreter()
    params = {'base': Path(os.getcwd()), 'Kmin': 2.5, 'Kmax': 12.5, 'kweight': 2.0, 'deltak': 0.05, 'rbkg': 1.0,
              'bkgkw': 2.0, 'bkgkmax': 15.0, 'front': ["tests/cu_test_files/cu_paths/path_75/feff"],
              'CSV': "tests/cu_test_files/cu_paths/cu_10k.xmu", 'optimize': False, 'series_index': 0, 'series': False}
    dirs = 'tests/cu_test_files/cu_results/'

    paths = [1, 2, 3, 5, 6, 10, 14, 16, 28, 30, 36, 40, 42]
    larch_score.larch_init(params['CSV'], params)

    Test_Result = EXAFS_Analysis.EXAFS_Analysis(paths, dirs, params)
    Test_Result.larch_init()

    def test_larch_extract_data_shape(self):
        """
        Test larch extract to have the correct matrix shape
        """
        TestCase.Test_Result.extract_data(verbose=False)
        result_shape = TestCase.Test_Result.bestFit_mat.shape
        self.assertEqual(result_shape[0], 13)
        self.assertEqual(result_shape[1], 4)

    def test_larch_extract_data_info(self):
        """
        Test larch extract to have the right data
        """
        TestCase.Test_Result.extract_data(verbose=False)
        test_result = TestCase.Test_Result.bestFit_mat
        real_result = np.array([[0.76, 0.03, 0.002, -0.02],
                                [0.22, 0.03, 0.005, 0.01],
                                [0.38, 0.03, 0.013, -0.05],
                                [0.86, 0.03, 0.001, 0.],
                                [0.54, 0.03, 0.008, 0.01],
                                [0.62, 0.03, 0.009, -0.05],
                                [0.94, 0.03, 0.013, 0.09],
                                [0.90, 0.03, 0.012, 0.09],
                                [0.82, 0.03, 0.001, 0.07],
                                [0.71, 0.03, 0.006, 0.06],
                                [0.53, 0.03, 0.011, -0.08],
                                [0.06, 0.03, 0.008, 0.10],
                                [0.50, 0.03, 0.006, 0.05]])

        self.assertTrue(np.allclose(test_result, real_result))

    def test_chi(self):
        """
        Test larch score data
        """
        TestCase.Test_Result.extract_data(verbose=False)
        TestCase.Test_Result.larch_score(verbose=False)

        Chi2 = np.round(TestCase.Test_Result.loss, 5)
        Chir2 = TestCase.Test_Result.chir2

        self.assertAlmostEqual(Chi2, 707.70957, 3)
        self.assertAlmostEqual(Chir2, 4.395712, 3)

    def test_bestFit_r(self):
        """Test best fit R data value
        """
        TestCase.Test_Result.extract_data(verbose=False)
        TestCase.Test_Result.larch_init()
        TestCase.Test_Result.larch_score(verbose=False)

        best_Fit_r = TestCase.Test_Result.best_Fit_r

        self.assertEqual(len(best_Fit_r), 13)
        best_Fit_r_1 = best_Fit_r[0]
        # S02
        self.assertAlmostEqual(best_Fit_r_1[0], 0.76, 3)
        self.assertAlmostEqual(best_Fit_r_1[1], 0.03,3)
        self.assertAlmostEqual(best_Fit_r_1[2], 0.002,3)
        # Kmin/Kmax/kweight - fixed by `params`, not by the GA fit.
        self.assertAlmostEqual(best_Fit_r_1[3], 2.5327,4)
        self.assertAlmostEqual(best_Fit_r_1[4], 12.0,3)
        self.assertAlmostEqual(best_Fit_r_1[5], 2.0,3)
        self.assertIsInstance(best_Fit_r_1[6], list)

    def test_larch_init(self):
        """
        Test individual export_paths shape
        """
        TestCase.Test_Result.extract_data(verbose=False)
        TestCase.Test_Result.larch_score(verbose=False)
        TestCase.Test_Result.individual_fit()

        ind_export_paths = TestCase.Test_Result.ind_export_paths
        self.assertEqual(ind_export_paths.shape[0], 26)
        self.assertEqual(ind_export_paths.shape[1], 401)
