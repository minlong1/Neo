import unittest

from PhysicsModules.EXAFS.exafs_neo.exafs_pop import NeoPopulations
from PhysicsModules.EXAFS.exafs_neo.neoPars import NeoPars
from PhysicsModules.EXAFS.exafs_neo.neoSolver import NeoSolver


def create_neosolver_operator(sol_opt):
    inputs_pars = {'data_file': 'tests/cu_test_files/cu_paths/cu_10k.xmu', 'output_file': '',
                   'feff_file': 'tests/cu_test_files/cu_paths/path_75/feff', 'kmin': 0.95,
                   'kmax': 9.775,
                   'kweight': 3.0, 'pathrange': [1, 2, 3, 4, 5],
                   'deltak': 0.05, 'rbkg': 1.1, 'bkgkw': 1.0, 'bkgkmax': 15.0,
                   'solOpt': sol_opt}

    exafs_Pars = NeoPars()
    exafs_Pars.read_inputs(inputs_pars)
    neo_population = NeoPopulations(exafs_Pars)
    neo_population.initialize_populations()

    return exafs_Pars


class TestNeoSolverBase(unittest.TestCase):
    inputs_pars = {'data_file': '../path_files/Cu/cu_10k.xmu', 'output_file': '',
                   'feff_file': '../path_files/Cu/path_75/feff', 'kmin': 0.95,
                   'kmax': 9.775,
                   'kweight': 3.0, 'pathrange': [1, 2, 3, 4, 5],
                   'deltak': 0.05, 'rbkg': 1.1, 'bkgkw': 1.0, 'bkgkmax': 15.0,
                   'mut_options': 1,
                   'solOpt': 0}

    def test_neosolver_base(self):
        pass


class TestNeoSolver(unittest.TestCase):

    def test_neosolver_ga(self):
        exafs_Pars = create_neosolver_operator(0)
        solver = NeoSolver()
        solver.initialize(exafs_pars=exafs_Pars)

        self.assertEqual(solver.solver_type, 0)
        self.assertEqual(solver.solver_operator.solver_type, 0)
        self.assertEqual(solver.solver_operator.solver_operator, 'Genetic Algorithm')

    def test_neosolver_ga_rech(self):
        exafs_Pars = create_neosolver_operator(1)
        solver = NeoSolver()
        solver.initialize(exafs_pars=exafs_Pars)

        self.assertEqual(solver.solver_type, 1)
        self.assertEqual(solver.solver_operator.solver_type, 1)
        self.assertEqual(solver.solver_operator.solver_operator, 'Genetic Algorithm with Rechenberg')

    def test_neosolver_de(self):
        exafs_Pars = create_neosolver_operator(2)
        solver = NeoSolver()
        solver.initialize(exafs_pars=exafs_Pars)

        self.assertEqual(solver.solver_type, 2)
        self.assertEqual(solver.solver_operator.solver_type, 2)
        self.assertEqual(solver.solver_operator.solver_operator, 'Differential Evolution')

    def test_neosolver_de_improves_fitness(self):
        import numpy as np

        np.random.seed(11)
        inputs_pars = {'data_file': 'tests/cu_test_files/cu_paths/cu_10k.xmu', 'output_file': '',
                       'feff_file': 'tests/cu_test_files/cu_paths/path_75/feff', 'kmin': 0.95,
                       'kmax': 9.775, 'kweight': 3.0, 'pathrange': [1, 2, 3, 4, 5],
                       'deltak': 0.05, 'rbkg': 1.1, 'bkgkw': 1.0, 'bkgkmax': 15.0,
                       'nPops': 20, 'solOpt': 2}
        exafs_Pars = NeoPars()
        exafs_Pars.read_inputs(inputs_pars)
        neo_population = NeoPopulations(exafs_Pars)
        neo_population.initialize_populations()
        before = exafs_Pars.bestFitPars.globBestVal

        solver = NeoSolver()
        solver.initialize(exafs_pars=exafs_Pars)
        for _ in range(5):
            solver.solve(neo_population, None, None, None, exafs_Pars)

        after = exafs_Pars.bestFitPars.globBestVal
        self.assertLess(after, before)


if __name__ == '__main__':
    unittest.main()
