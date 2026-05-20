"""
Unit test for EXAFS Analysis
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
        real_result = np.array([[0.95, 0.13, 0.004, -0.02],
                                [0.92, 0.13, 0.006, -0.03],
                                [0.55, 0.13, 0.002, -0.],
                                [0.84, 0.13, 0.004, -0.01],
                                [0.88, 0.13, 0.005, 0.06],
                                [0.55, 0.13, 0.005, 0.04],
                                [0.95, 0.13, 0.011, -0.06],
                                [0.35, 0.13, 0.001, 0.01],
                                [0.22, 0.13, 0.008, -0.07],
                                [0.95, 0.13, 0.006, -0.09],
                                [0.66, 0.13, 0.001, -0.06],
                                [0.75, 0.13, 0.001, 0.01],
                                [0.85, 0.13, 0.004, -0.08]])

        self.assertTrue(np.allclose(test_result, real_result))

    def test_chi(self):
        """
        Test larch score data
        """
        TestCase.Test_Result.extract_data(verbose=False)
        TestCase.Test_Result.larch_score(verbose=False)

        Chi2 = np.round(TestCase.Test_Result.loss, 5)
        Chir2 = TestCase.Test_Result.chir2

        self.assertAlmostEqual(Chi2, 52.562, 3)
        self.assertAlmostEqual(Chir2, 0.326471, 3)

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
        self.assertAlmostEqual(best_Fit_r_1[0], 0.95, 3)
        self.assertAlmostEqual(best_Fit_r_1[1], 0.13,3)
        self.assertAlmostEqual(best_Fit_r_1[2], 0.004,3)
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
