import unittest

import numpy as np

from Solvers.core import (
    GeneRange,
    Individual,
    ParameterSpace,
    Population,
    RunState,
    SolverResult,
    local_search_refine,
)
from Solvers.tests.toy_problem import QuadraticProblem


class TestGeneRange(unittest.TestCase):
    def test_bounds(self):
        gene = GeneRange(np.arange(0.0, 1.01, 0.1))
        self.assertAlmostEqual(gene.low, 0.0)
        self.assertAlmostEqual(gene.high, 1.0)

    def test_sample_within_values(self):
        gene = GeneRange([1.0, 2.0, 3.0])
        for _ in range(20):
            self.assertIn(gene.sample(), [1.0, 2.0, 3.0])

    def test_clip(self):
        gene = GeneRange([0.0, 0.5, 1.0])
        self.assertEqual(gene.clip(2.0), 1.0)
        self.assertEqual(gene.clip(-1.0), 0.0)

    def test_empty_rejected(self):
        with self.assertRaises(ValueError):
            GeneRange([])

    def test_from_bounds(self):
        gene = GeneRange.from_bounds(0.0, 1.0, 0.25)
        self.assertEqual(len(gene), 4)


class TestParameterSpace(unittest.TestCase):
    def setUp(self):
        self.space = ParameterSpace(
            [GeneRange([0.0, 1.0]), GeneRange([-1.0, 0.0, 1.0]), GeneRange([5.0])]
        )

    def test_n_genes(self):
        self.assertEqual(self.space.n_genes, 3)
        self.assertEqual(len(self.space), 3)

    def test_sample_shape_and_membership(self):
        genes = self.space.sample()
        self.assertEqual(genes.shape, (3,))
        self.assertIn(genes[0], [0.0, 1.0])
        self.assertEqual(genes[2], 5.0)

    def test_clip(self):
        clipped = self.space.clip(np.array([9.0, -9.0, 5.0]))
        self.assertEqual(list(clipped), [1.0, -1.0, 5.0])

    def test_limits(self):
        self.assertEqual(self.space.limits(1), (-1.0, 1.0))


class TestIndividual(unittest.TestCase):
    def setUp(self):
        self.space = ParameterSpace([GeneRange(np.arange(0, 1.01, 0.1)) for _ in range(5)])

    def test_random_init(self):
        ind = Individual(self.space)
        self.assertEqual(len(ind), 5)

    def test_explicit_genes(self):
        ind = Individual(self.space, np.array([0.1, 0.2, 0.3, 0.4, 0.5]))
        self.assertAlmostEqual(ind.get_gene(2), 0.3)

    def test_wrong_shape_rejected(self):
        with self.assertRaises(ValueError):
            Individual(self.space, np.zeros(3))

    def test_mutate_always(self):
        ind = Individual(self.space)
        n = ind.mutate(chance=1.0)
        self.assertEqual(n, 5)

    def test_mutate_never(self):
        ind = Individual(self.space)
        before = ind.get()
        n = ind.mutate(chance=0.0)
        self.assertEqual(n, 0)
        np.testing.assert_array_equal(ind.get(), before)

    def test_copy_is_independent(self):
        ind = Individual(self.space)
        clone = ind.copy()
        clone.set_gene(0, 99.0)
        self.assertNotEqual(ind.get_gene(0), 99.0)


class TestPopulation(unittest.TestCase):
    def setUp(self):
        self.problem = QuadraticProblem()
        self.state = RunState(n_gen=10)
        self.pop = Population(self.problem, self.state, n_pops=20)
        self.pop.initialize_populations()

    def test_initialize(self):
        self.assertEqual(len(self.pop), 20)
        self.assertEqual(len(self.pop.population_sorted), 20)

    def test_sorted_ascending(self):
        scores = [s for _, s in self.pop.population_sorted]
        self.assertEqual(scores, sorted(scores))

    def test_best_tracked_in_state(self):
        self.assertIsNotNone(self.state.globBestInd)
        self.assertEqual(self.state.currBestVal, self.pop.population_sorted[0][1])
        self.assertLessEqual(self.state.globBestVal, self.state.currBestVal)


class TestRunState(unittest.TestCase):
    def test_update_best_improvement(self):
        state = RunState(n_gen=10)
        state.update_best("ind1", 5.0)
        self.assertEqual(state.globBestVal, 5.0)
        state.update_best("ind2", 3.0)
        self.assertEqual(state.globBestInd, "ind2")
        self.assertEqual(state.globBestVal, 3.0)

    def test_update_best_no_regression(self):
        state = RunState(n_gen=10)
        state.update_best("ind1", 3.0)
        state.update_best("ind2", 5.0)
        self.assertEqual(state.globBestInd, "ind1")
        self.assertEqual(state.globBestVal, 3.0)

    def test_small_delta_zeroes_bestdiff(self):
        state = RunState(n_gen=10)
        state.update_best("a", 3.0)
        state.update_best("b", 3.005)
        self.assertEqual(state.bestDiff, 0)

    def test_second_half(self):
        state = RunState(n_gen=10)
        self.assertFalse(state.second_half)
        state.currGen = 5
        self.assertTrue(state.second_half)


class TestLocalSearchRefine(unittest.TestCase):
    """Solvers.core.refinement.local_search_refine — the genome-agnostic
    generalization of Solvers.demcmc.refinement.local_search_refine used
    by PSOSolver's optional local_search option."""

    def setUp(self):
        self.problem = QuadraticProblem()  # target in [-0.5, 0.5], genes in [-1, 1]
        self.state = RunState(n_gen=10)
        self.pop = Population(self.problem, self.state, n_pops=10)
        self.pop.initialize_populations()

    def test_rejects_when_population_min_unbeatable(self):
        np.random.seed(0)
        before = [ind.genes.copy() for ind in self.pop.population]
        n_updated = local_search_refine(self.pop, self.problem, population_min=-1e9, n_samples=20)
        self.assertEqual(n_updated, 0)
        for ind, snapshot in zip(self.pop.population, before):
            np.testing.assert_array_equal(ind.genes, snapshot)

    def test_accepts_genuine_improvement(self):
        np.random.seed(1)
        bad_genes = np.array([1.0, 1.0, 1.0, 1.0])  # far from target
        individual = Individual(self.problem.space, bad_genes)
        bad_score = self.problem.fitness(bad_genes)

        n_updated = local_search_refine(
            self.pop,
            self.problem,
            population_min=bad_score,
            individuals=[individual],
            n_samples=200,
            perturb_scale=0.5,
        )
        self.assertEqual(n_updated, 1)
        self.assertLess(self.problem.fitness(individual.genes), bad_score)

    def test_respects_bounds(self):
        np.random.seed(2)
        edge_genes = np.array([1.0, 1.0, 1.0, 1.0])  # already at the upper bound
        individual = Individual(self.problem.space, edge_genes)

        local_search_refine(
            self.pop,
            self.problem,
            population_min=1e9,  # trivially beatable
            individuals=[individual],
            n_samples=50,
            perturb_scale=1.0,
        )
        for i in range(len(individual.genes)):
            low, high = self.problem.space.limits(i)
            self.assertGreaterEqual(individual.genes[i], low)
            self.assertLessEqual(individual.genes[i], high)


class TestSolverResult(unittest.TestCase):
    def test_collect_and_roundtrip(self):
        import os
        import tempfile

        state = RunState(n_gen=2)
        state.update_best("best", 1.5)
        result = SolverResult()
        result.collect(state, None)
        self.assertEqual(result.historyBest, [1.5])

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "result.pkl")
            result.save(path)
            loaded = SolverResult.load(path)
            self.assertEqual(loaded.historyBest, [1.5])
            self.assertEqual(loaded.best_value, 1.5)


if __name__ == "__main__":
    unittest.main()
