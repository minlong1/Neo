import unittest

from PhysicsModules.EXAFS.exafs_neo.exafs_pop import NeoPopulations
from PhysicsModules.EXAFS.exafs_neo.neoPars import NeoPars
from PhysicsModules.EXAFS.exafs_neo.neoCrossOver import NeoCrossover


def create_test_crossover_operator(croOpt):
    inputs_pars = {'data_file': 'tests/cu_test_files/cu_paths/cu_10k.xmu', 'output_file': '',
                   'feff_file': 'tests/cu_test_files/cu_paths/path_75/feff', 'kmin': 0.95,
                   'kmax': 9.775,
                   'npop': 2,
                   'kweight': 3.0, 'pathrange': [1, 2, 3, 4, 5],
                   'deltak': 0.05, 'rbkg': 1.1, 'bkgkw': 1.0, 'bkgkmax': 15.0,
                   'mut_options': 1,
                   'croOpt': croOpt}
    exafs_Pars = NeoPars()
    exafs_Pars.read_inputs(inputs_pars)
    neo_population = NeoPopulations(exafs_Pars)
    neo_population.initialize_populations()

    crossover_operator = NeoCrossover()
    crossover_operator.initialize(exafs_pars=exafs_Pars)

    return crossover_operator, neo_population, (neo_population[0][0], neo_population[1][0])


class TestNeoCrossOverBase(unittest.TestCase):
    inputs_pars = {'data_file': '../path_files/Cu/cu_10k.xmu', 'output_file': '',
                   'feff_file': '../path_files/Cu/path_75/feff', 'kmin': 0.95,
                   'kmax': 9.775,
                   'kweight': 3.0, 'pathrange': [1, 2, 3, 4, 5],
                   'deltak': 0.05, 'rbkg': 1.1, 'bkgkw': 1.0, 'bkgkmax': 15.0,
                   'mut_options': 1,
                   'croOpt': 0}

    def test_neocrossover_base(self):
        pass


class TestNeoCrossOver(unittest.TestCase):
    def test_neocrossover(self):
        crossover_operator, neo_population, (ind1, ind2) = create_test_crossover_operator(0)

        self.assertEqual(crossover_operator.crossover_type, 0)
        self.assertEqual(crossover_operator.crossover_operator.croOpt, 0)
        self.assertEqual(crossover_operator.crossover_operator.croType, 'Uniform Crossover')


class TestNeoCrossOverFunc(unittest.TestCase):

    def test_neocrossover_uniform(self):
        num_opt = 0
        crossover_operator, neo_population, (ind1, ind2) = create_test_crossover_operator(num_opt)

        self.assertEqual(crossover_operator.crossover_type, num_opt)
        self.assertEqual(crossover_operator.crossover_operator.croOpt, num_opt)
        self.assertEqual(crossover_operator.crossover_operator.croType, 'Uniform Crossover')

        new_ind = crossover_operator.crossover_single(neo_population, ind1, ind2)
        for path in range(5):
            for par in range(4):
                self.assertTrue(
                    ind1.get()[path][par] == new_ind.get()[path][par] or ind2.get()[path][par] == new_ind.get()[path][
                        par])

    def test_neocrossover_single_point(self):
        num_opt = 1
        crossover_operator, neo_population, (ind1, ind2) = create_test_crossover_operator(num_opt)

        self.assertEqual(crossover_operator.crossover_type, num_opt)
        self.assertEqual(crossover_operator.crossover_operator.croOpt, num_opt)
        self.assertEqual(crossover_operator.crossover_operator.croType, 'Single Point Crossover')
        # print(ind1)
        # print(ind2)

        new_ind = crossover_operator.crossover_single(neo_population, ind1, ind2)
        # print(new_ind)

        for path in range(5):
            for par in range(4):
                self.assertTrue(
                    ind1.get()[path][par] == new_ind.get()[path][par] or ind2.get()[path][par] == new_ind.get()[path][
                        par])

    def test_neocrossover_arithmetic_crossover(self):
        num_opt = 3
        crossover_operator, neo_population, (ind1, ind2) = create_test_crossover_operator(num_opt)

        self.assertEqual(crossover_operator.crossover_type, num_opt)
        self.assertEqual(crossover_operator.crossover_operator.croOpt, num_opt)
        self.assertEqual(crossover_operator.crossover_operator.croType, 'Arithmetic Crossover')

        new_ind = crossover_operator.crossover_single(neo_population, ind1, ind2)
        for path in range(5):
            for par in range(4):
                self.assertTrue(
                    ind1.get()[path][par] == new_ind.get()[path][par] or ind2.get()[path][par] == new_ind.get()[path][
                        par])

    def test_neocrossover_or_crossover(self):
        num_opt = 4
        crossover_operator, neo_population, (ind1, ind2) = create_test_crossover_operator(num_opt)

        self.assertEqual(crossover_operator.crossover_type, num_opt)
        self.assertEqual(crossover_operator.crossover_operator.croOpt, num_opt)
        self.assertEqual(crossover_operator.crossover_operator.croType, 'Or Crossover')

        new_ind = crossover_operator.crossover_single(neo_population, ind1, ind2)
        for path in range(5):
            for par in range(4):
                self.assertTrue(
                    ind1.get()[path][par] == new_ind.get()[path][par] or ind2.get()[path][par] == new_ind.get()[path][
                        par])

    def test_neocrossover_average_crossover(self):
        num_opt = 5
        crossover_operator, neo_population, (ind1, ind2) = create_test_crossover_operator(num_opt)

        self.assertEqual(crossover_operator.crossover_type, num_opt)
        self.assertEqual(crossover_operator.crossover_operator.croOpt, num_opt)
        self.assertEqual(crossover_operator.crossover_operator.croType, 'Average Crossover')

        new_ind = crossover_operator.crossover_single(neo_population, ind1, ind2)

        for path in range(5):
            for par in [0, 2, 3]:
                average = (ind1.get()[path][par] + ind2.get()[path][par]) / 2
                self.assertEqual(average, new_ind.get()[path][par])


if __name__ == '__main__':
    unittest.main()
