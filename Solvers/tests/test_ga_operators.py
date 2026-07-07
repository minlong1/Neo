import unittest

import numpy as np

from Solvers.core import Population, RunState
from Solvers.ga import (
    AverageCrossover,
    BoundedMutator,
    MetropolisMutator,
    PerGeneMutator,
    PerIndividualMutator,
    RouletteWheelSelector,
    SinglePointCrossover,
    UniformCrossover,
    create_crossover,
    create_mutator,
    create_selector,
)
from Solvers.tests.toy_problem import QuadraticProblem


def make_population(n_pops=20, n_gen=10):
    problem = QuadraticProblem()
    state = RunState(n_gen=n_gen)
    state.currGen = 1
    pop = Population(problem, state, n_pops=n_pops)
    pop.initialize_populations()
    return pop


class TestSelectors(unittest.TestCase):
    def test_roulette_wheel_counts(self):
        pop = make_population(n_pops=20)
        selector = RouletteWheelSelector(n_pops=20, n_best_sample=0.3, n_lucky_sample=0.2)
        selector.select(pop)
        # 6 best + 4 lucky
        self.assertEqual(len(pop.next_population), 10)
        self.assertEqual(selector.nCross, 10)

    def test_best_survive(self):
        pop = make_population(n_pops=20)
        best = pop.population_sorted[0][0]
        selector = RouletteWheelSelector(n_pops=20)
        selector.select(pop)
        self.assertIn(best, pop.next_population)

    def test_factory_invalid(self):
        with self.assertRaises(ValueError):
            create_selector(99, n_pops=10)


class TestCrossovers(unittest.TestCase):
    def setUp(self):
        self.pop = make_population(n_pops=20)
        self.ind1 = self.pop.population[0]
        self.ind2 = self.pop.population[1]

    def assert_child_mixes_parents(self, child):
        for i in range(len(child)):
            self.assertIn(
                child.genes[i], [self.ind1.genes[i], self.ind2.genes[i]]
            )

    def test_uniform(self):
        child = UniformCrossover(n_cross=1).crossover_pair(self.pop, self.ind1, self.ind2)
        self.assert_child_mixes_parents(child)

    def test_single_point(self):
        child = SinglePointCrossover(n_cross=1).crossover_pair(
            self.pop, self.ind1, self.ind2, co_point=2
        )
        np.testing.assert_array_equal(child.genes[:2], self.ind1.genes[:2])
        np.testing.assert_array_equal(child.genes[2:], self.ind2.genes[2:])

    def test_average(self):
        child = AverageCrossover(n_cross=1).crossover_pair(self.pop, self.ind1, self.ind2)
        np.testing.assert_allclose(child.genes, (self.ind1.genes + self.ind2.genes) / 2)

    def test_crossover_refills_population(self):
        selector = RouletteWheelSelector(n_pops=20)
        selector.select(self.pop)
        crossover = UniformCrossover(n_cross=selector.nCross)
        crossover.crossover(self.pop)
        self.assertEqual(len(self.pop.population), 20)

    def test_factory_invalid(self):
        with self.assertRaises(ValueError):
            create_crossover(99, n_cross=1)


class TestMutators(unittest.TestCase):
    def test_per_individual_all(self):
        pop = make_population(n_pops=10)
        originals = list(pop.population)
        PerIndividualMutator(mut_chance=1.0).mutate(pop)
        self.assertTrue(all(ind not in originals for ind in pop.population))

    def test_per_gene_none(self):
        pop = make_population(n_pops=10)
        before = [ind.get() for ind in pop.population]
        PerGeneMutator(mut_chance=0.0).mutate(pop)
        for ind, genes in zip(pop.population, before):
            np.testing.assert_array_equal(ind.get(), genes)

    def test_metropolis_runs(self):
        pop = make_population(n_pops=10)
        pop.state.currGen = 2
        pop.state.bestDiff = 1.0
        MetropolisMutator(mut_chance=1.0).mutate(pop)
        self.assertEqual(len(pop.population), 10)

    def test_bounded_stays_in_bounds(self):
        pop = make_population(n_pops=10)
        BoundedMutator(mut_chance=1.0).mutate(pop)
        for ind in pop.population:
            for j in range(len(ind)):
                low, high = pop.problem.space.limits(j)
                self.assertGreaterEqual(ind.get_gene(j), low - 1e-9)
                self.assertLessEqual(ind.get_gene(j), high + 1e-9)

    def test_factory_invalid(self):
        with self.assertRaises(ValueError):
            create_mutator(99)


if __name__ == "__main__":
    unittest.main()
