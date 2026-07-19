import unittest

from PhysicsModules.EXAFS.exafs_neo.exafs_pop import NeoPopulations
from PhysicsModules.EXAFS.exafs_neo.neoPars import NeoPars
from PhysicsModules.EXAFS.exafs_neo.neoMutator import NeoMutator


def create_neo_mutator_operator(mut_opt):
    inputs_pars = {'data_file': 'tests/cu_test_files/cu_paths/cu_10k.xmu', 'output_file': '',
                   'feff_file': 'tests/cu_test_files/cu_paths/path_75/feff', 'kmin': 0.95,
                   'kmax': 9.775,
                   'kweight': 3.0, 'pathrange': [1, 2, 3, 4, 5],
                   'deltak': 0.05, 'rbkg': 1.1, 'bkgkw': 1.0, 'bkgkmax': 15.0,
                   'mut_options': mut_opt,
                   'croOpt': 0}
    exafs_Pars = NeoPars()
    exafs_Pars.read_inputs(inputs_pars)
    neo_population = NeoPopulations(exafs_Pars)
    neo_population.initialize_populations()

    return exafs_Pars


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


class TestNeoMutator(unittest.TestCase):

    def test_neomutator_per_individual(self):
        exafs_Pars = create_neo_mutator_operator(0)
        mutator_operator = NeoMutator()
        mutator_operator.initialize(exafs_pars=exafs_Pars)

        self.assertEqual(mutator_operator.mutator_type, 0)
        self.assertEqual(mutator_operator.mutator.mutOpt, 1)
        self.assertEqual(mutator_operator.mutator.mutType, 'Mutate Per Individual')

    def test_neomutator_per_path(self):
        exafs_Pars = create_neo_mutator_operator(1)
        mutator_operator = NeoMutator()
        mutator_operator.initialize(exafs_pars=exafs_Pars)

        self.assertEqual(mutator_operator.mutator_type, 1)
        self.assertEqual(mutator_operator.mutator.mutOpt, 2)
        self.assertEqual(mutator_operator.mutator.mutType, 'Mutate Per Path')

    def test_neomutator_per_trait(self):
        exafs_Pars = create_neo_mutator_operator(2)
        mutator_operator = NeoMutator()
        mutator_operator.initialize(exafs_pars=exafs_Pars)

        self.assertEqual(mutator_operator.mutator_type, 2)
        self.assertEqual(mutator_operator.mutator.mutOpt, 3)
        self.assertEqual(mutator_operator.mutator.mutType, 'Mutate Per Trait')

    def test_neomutator_metropolis(self):
        exafs_Pars = create_neo_mutator_operator(3)
        mutator_operator = NeoMutator()
        mutator_operator.initialize(exafs_pars=exafs_Pars)

        self.assertEqual(mutator_operator.mutator_type, 3)
        self.assertEqual(mutator_operator.mutator.mutOpt, 4)
        self.assertEqual(mutator_operator.mutator.mutType, 'Mutate Metropolis')

    def test_neomutator_bounded(self):
        exafs_Pars = create_neo_mutator_operator(4)
        mutator_operator = NeoMutator()
        mutator_operator.initialize(exafs_pars=exafs_Pars)

        self.assertEqual(mutator_operator.mutator_type, 4)
        self.assertEqual(mutator_operator.mutator.mutOpt, 5)
        self.assertEqual(mutator_operator.mutator.mutType, 'Mutate Bounded Per Range')

    def test_neomutator_de(self):
        exafs_Pars = create_neo_mutator_operator(5)
        mutator_operator = NeoMutator()
        mutator_operator.initialize(exafs_pars=exafs_Pars)

        self.assertEqual(mutator_operator.mutator_type, 5)
        self.assertEqual(mutator_operator.mutator.mutOpt, 6)
        self.assertEqual(mutator_operator.mutator.mutType, 'Mutate DE')

    def test_neomutator_resyncs_mutchance_every_call(self):
        """mutate() must re-read mutPars.mutChance each call, not just at
        initialize() — this is what lets GA_Rechenberg's adaptive update
        (which writes mutPars.mutChance) actually reach the operator."""
        exafs_Pars = create_neo_mutator_operator(0)
        mutator_operator = NeoMutator()
        mutator_operator.initialize(exafs_pars=exafs_Pars)
        neo_population = NeoPopulations(exafs_Pars)
        neo_population.initialize_populations()

        self.assertEqual(mutator_operator.mutator.mut_chance, exafs_Pars.mutPars.mutChance)

        exafs_Pars.mutPars.mutChance = 0.85
        mutator_operator.mutate(neo_population)

        self.assertEqual(mutator_operator.mutator.mut_chance, 0.85)
        self.assertEqual(mutator_operator.mutator.mutChance, 0.85)


if __name__ == '__main__':
    unittest.main()
