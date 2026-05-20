import unittest

from PhysicsModules.EXAFS.exafs_neo import pathrange


class Test_Pathrange(unittest.TestCase):
    temp_obj = pathrange.Pathrange_limits([1, 2, 3, 4])

    def test_get_lim(self):
        temp_range = Test_Pathrange.temp_obj.get_lim()
        self.assertEqual(temp_range[0], [1, 2, 3, 4])
        self.assertAlmostEqual(temp_range[1][0], 0.0)
        self.assertAlmostEqual(temp_range[1][1], 0.95)
        self.assertAlmostEqual(temp_range[1][2], 0.01)

        self.assertAlmostEqual(temp_range[2][0], 0.001)
        self.assertAlmostEqual(temp_range[2][1], 0.015)
        self.assertAlmostEqual(temp_range[2][2], 0.001)

        self.assertAlmostEqual(temp_range[3][0], -0.1)
        self.assertAlmostEqual(temp_range[3][1], 0.1)
        self.assertAlmostEqual(temp_range[3][2], 0.01)
        # self.assertSetEqual()

    def test_get_lim_S02(self):
        temp_S02 = Test_Pathrange.temp_obj.get_lim_S02()
        self.assertAlmostEqual(temp_S02[0], 0.0)
        self.assertAlmostEqual(temp_S02[1], 0.95)
        self.assertAlmostEqual(temp_S02[2], 0.01)

    def test_get_lim_Sigma2(self):
        temp_Sigma2 = Test_Pathrange.temp_obj.get_lim_Sigma2()
        self.assertAlmostEqual(temp_Sigma2[0], 0.001)
        self.assertAlmostEqual(temp_Sigma2[1], 0.015)
        self.assertAlmostEqual(temp_Sigma2[2], 0.001)

    def test_get_lim_DeltaR(self):
        temp_deltaR = Test_Pathrange.temp_obj.get_lim_DeltaR()
        self.assertAlmostEqual(temp_deltaR[0], -0.1)
        self.assertAlmostEqual(temp_deltaR[1], 0.1)
        self.assertAlmostEqual(temp_deltaR[2], 0.01)
