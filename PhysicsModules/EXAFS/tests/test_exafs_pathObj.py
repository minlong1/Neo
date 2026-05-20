import unittest

import numpy as np

from PhysicsModules.EXAFS.exafs_neo import pathObj, pathrange


class Test_PathObj(unittest.TestCase):
    e0 = 0.37
    pathrange_Obj = pathrange.Pathrange_limits(0)
    pathObj = pathObj.PathObject(pathrange_Obj, 0.37)

    def test_get_shape(self):
        """Test shape of the get method
        """
        var_result = np.array(Test_PathObj.pathObj.get())
        self.assertEqual(var_result.shape[0], 4)

    def test_get_val(self):
        """Test the value of get
        """
        var_result = np.array(Test_PathObj.pathObj.get())
        self.assertEqual(var_result[1], 0.37)

    def test_getvar_shape(self):
        """Test shape of the get_var method
        """
        var_result = np.array(Test_PathObj.pathObj.get_var())
        self.assertEqual(var_result.shape[0], 3)


    def test_set_val(self):
        """Test the set method
        """
        pars = Test_PathObj.pathObj.get()
        Test_PathObj.pathObj.set(1, 2, 3)
        self.assertEqual(Test_PathObj.pathObj.s02, 1)
        self.assertEqual(Test_PathObj.pathObj.sigma2, 2)
        self.assertEqual(Test_PathObj.pathObj.deltaR, 3)
        Test_PathObj.pathObj.set(pars[0], pars[1], pars[2])

    def test_set_s02(self):
        """Test the set_s02 method
        """
        oldS02 = Test_PathObj.pathObj.s02
        Test_PathObj.pathObj.set_s02(1)
        self.assertEqual(Test_PathObj.pathObj.s02, 1)
        Test_PathObj.pathObj.set_s02(oldS02)

    def test_set_e0(self):
        """Test the set_e0 method
        """
        oldE0 = Test_PathObj.pathObj.e0
        Test_PathObj.pathObj.set_e0(1)
        self.assertEqual(Test_PathObj.pathObj.e0, 1)
        Test_PathObj.pathObj.set_e0(oldE0)

    def test_set_sigma2(self):
        """Test the set_sigma2 method
        """
        oldSigma2 = Test_PathObj.pathObj.sigma2
        Test_PathObj.pathObj.set_sigma2(1)
        self.assertEqual(Test_PathObj.pathObj.sigma2, 1)
        Test_PathObj.pathObj.set_sigma2(oldSigma2)

    def test_set_deltaR(self):
        """Test the set_deltaR method
        """
        oldDeltaR = Test_PathObj.pathObj.deltaR
        Test_PathObj.pathObj.set_deltaR(1)
        self.assertEqual(Test_PathObj.pathObj.deltaR, 1)
        Test_PathObj.pathObj.set_deltaR(oldDeltaR)
