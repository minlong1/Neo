import unittest

import numpy as np

from PhysicsModules.EXAFS.exafs_neo.neoPars import NeoPars, NeoMutPars


class TestEXAFSpars(unittest.TestCase):
    inputs_pars = {'data_file': 'tests/cu_test_files/cu_paths/cu_10k.xmu',
                   'output_file': '',
                   'feff_file': 'tests/cu_test_files/cu_paths/path_75/feff',
                   'kmin': 0.95,
                   'kmax': 9.775,
                   'kweight': 3.0, 'pathrange': [1, 2, 3, 4, 5],
                   'deltak': 0.05, 'rbkg': 1.2, 'bkgkw': 1.0, 'bkgkmax': 15.0}
    exafs_NeoPars = NeoPars()

    exafs_NeoPars.read_inputs(inputs_pars)

    def test_EXAFS_pars(self):
        self.assertEqual(TestEXAFSpars.exafs_NeoPars.exafsPars.kmin, 0.95)
        self.assertEqual(TestEXAFSpars.exafs_NeoPars.exafsPars.kmax, 9.775)
        self.assertEqual(TestEXAFSpars.exafs_NeoPars.exafsPars.dk, 0.05)
        self.assertEqual(TestEXAFSpars.exafs_NeoPars.exafsPars.kweight, 3.0)
        self.assertEqual(TestEXAFSpars.exafs_NeoPars.exafsPars.rbkg, 1.2)

    def test_EXAFS_notEq(self):
        self.assertNotEqual(TestEXAFSpars.exafs_NeoPars.exafsPars.kmin, 0.96)
        self.assertNotEqual(TestEXAFSpars.exafs_NeoPars.exafsPars.kmax, 9.776)
        self.assertNotEqual(TestEXAFSpars.exafs_NeoPars.exafsPars.dk, 0.06)
        self.assertNotEqual(TestEXAFSpars.exafs_NeoPars.exafsPars.kweight, 3.1)
        self.assertNotEqual(TestEXAFSpars.exafs_NeoPars.exafsPars.rbkg, 1.3)

    def test_EXAFS_calc_pars(self):
        self.assertEqual(TestEXAFSpars.exafs_NeoPars.exafsPars.small, 18)
        self.assertEqual(TestEXAFSpars.exafs_NeoPars.exafsPars.big, 195)
        self.assertEqual(TestEXAFSpars.exafs_NeoPars.exafsPars.mid, 178)


class TestMutPars(unittest.TestCase):

    def test_mutPars(self):
        mutPars = NeoMutPars()
        pass


class TestNeoPars(unittest.TestCase):
    inputs_pars = {'data_file': 'tests/cu_test_files/cu_paths/cu_10k.xmu',
                   'output_file': 'test_output.csv',
                   'feff_file': 'tests/cu_test_files/cu_paths/path_75/feff',
                   'kmin': 0.95,
                   'kmax': 9.775,
                   'kweight': 3.0,
                   'pathrange': [1, 2, 3, 4, 5],
                   'deltak': 0.05, 'rbkg': 1.2, 'bkgkw': 1.0, 'bkgkmax': 15.0, 'nPops': 150, 'nGen': 200,
                   'mutChance': 0.15}

    exafs_NeoPars = NeoPars()
    exafs_NeoPars.read_inputs(inputs_pars)

    def test_runPars(self):
        self.assertEqual(TestEXAFSpars.exafs_NeoPars.runPars.currGen, 1)

    def test_fixedPars(self):
        self.assertEqual(TestNeoPars.exafs_NeoPars.fixedPars.nPops, 150)
        self.assertEqual(TestNeoPars.exafs_NeoPars.fixedPars.nGen, 200)
        self.assertEqual(TestNeoPars.exafs_NeoPars.fixedPars.selOpt, 0)

    def test_filePars(self):
        # self.assertEqual(TestNeoPars.exafs_NeoPars.neoFilePars.data_file, 'cu_test_files/cu_paths/cu_10k.xmu')
        self.assertEqual(TestNeoPars.exafs_NeoPars.neoFilePars.output_file, 'test_output.csv')
        # self.assertEqual(TestNeoPars.exafs_NeoPars.neoFilePars.feff_file, 'cu_test_files/cu_paths/path_75/feff')

    def test_mutPars(self):
        self.assertEqual(TestNeoPars.exafs_NeoPars.mutPars.mutOpt, 1)
        self.assertEqual(TestNeoPars.exafs_NeoPars.mutPars.mutChance, 0.15)
        self.assertEqual(TestNeoPars.exafs_NeoPars.mutPars.mutChanceE0, 0.3)

    def test_crossPars(self):
        self.assertEqual(TestNeoPars.exafs_NeoPars.crossPars.croOpt, 0)

    def test_selPars(self):
        self.assertEqual(TestNeoPars.exafs_NeoPars.selPars.selOpt, 0)

    def test_bestFitPars(self):
        self.assertEqual(TestNeoPars.exafs_NeoPars.bestFitPars.currBestVal, np.inf)
        self.assertEqual(TestNeoPars.exafs_NeoPars.bestFitPars.globBestVal, np.inf)
        self.assertEqual(TestNeoPars.exafs_NeoPars.bestFitPars.bestDiff, np.inf)
