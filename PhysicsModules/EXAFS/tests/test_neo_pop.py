import unittest
import numpy as np
from PhysicsModules.EXAFS.exafs_neo.exafs_pop import NeoPopulations
from PhysicsModules.EXAFS.exafs_neo.neoPars import NeoPars


class MyTestCase(unittest.TestCase):

    def test_neoPars(self):
        inputs_pars = {'data_file': 'tests/cu_test_files/cu_paths/cu_10k.xmu',
                       'output_file': '',
                       'feff_file': 'tests/cu_test_files/cu_paths/path_75/feff',
                       'kmin': 0.95,
                       'kmax': 9.775,
                       'kweight': 3.0, 'pathrange': [1, 2, 3, 4, 5],
                       'deltak': 0.05, 'rbkg': 1.1, 'bkgkw': 1.0, 'bkgkmax': 15.0}
        exafs_NeoPars = NeoPars()
        exafs_NeoPars.read_inputs(inputs_pars)
        neo_population = NeoPopulations(exafs_NeoPars)
        neo_population.initialize_populations()
        self.assertIsInstance(neo_population[0], tuple)
        # self.assertIsInstance(neo_population[0][0], exafs_test.individual.Individual)

    def test_neoPops(self):
        inputs_pars = {'data_file': 'tests/cu_test_files/cu_paths/cu_10k.xmu',
                       'output_file': '',
                       'feff_file': 'tests/cu_test_files/cu_paths/path_75/feff',
                       'kmin': 3.0,
                       'kmax': 17.0,
                       'kweight': 2.0,
                       'pathrange': [1, 2, 3, 5, 6, 10, 14, 16, 28, 30, 36, 40, 42],
                       'deltak': 0.05, 'rbkg': 1.2, 'bkgkw': 1.0, 'bkgkmax': 25.0,
                       'individual':True}
        # Your original list
        raw_data =  np.array([
            [0.94, 1.33, 0.004, -0.01],
            [0.78, 1.33, 0.005, -0.02],
            [0.5, 1.33, 0.007, -0.01],
            [0.95, 1.33, 0.005, 0.02],
            [0.37, 1.33, 0.003, -0.09],
            [0.27, 1.33, 0.004, 0.03],
            [0.69, 1.33, 0.006, -0.05],
            [0.45, 1.33, 0.005, -0.02],
            [0.7, 1.33, 0.01, 0.07],
            [0.82, 1.33, 0.009, 0.02],
            [0.62, 1.33, 0.006, -0.04],
            [0.35, 1.33, 0.005, -0.06],
            [0.81, 1.33, 0.004, -0.08]
        ])
        exafs_NeoPars = NeoPars()
        exafs_NeoPars.read_inputs(inputs_pars)
        neo_population = NeoPopulations(exafs_NeoPars)
        individual = neo_population.generate_individual()

        for i,path in enumerate(raw_data):
            individual.set_path(i,path[0],path[2],path[3])
        individual.set_e0(raw_data[0][1])

        neo_population.population.append(individual)
        self.assertAlmostEqual(neo_population.eval_population()[0],398.76,2)
        self.assertEqual(len(neo_population.population),1)





if __name__ == '__main__':
    unittest.main()
