import unittest

import numpy as np

from Solvers import DEMCMCSolver, get_solver
from Solvers.demcmc import BNNRegressionProblem, MLPStructure, PosteriorResult
from Solvers.demcmc.hyper_mutation import StagnationTracker, sample_hyperparameters
from Solvers.demcmc.mutation import (
    MUTATION_OPERATORS,
    add_mcmc_noise,
    mutate_best_1,
    mutate_best_2,
    mutate_rand_1,
    mutate_rand_2,
)
from Solvers.demcmc.refinement import local_search_refine, svd_refine


def make_regression_data(n=40, seed=0):
    rng = np.random.RandomState(seed)
    X = np.linspace(-3, 3, n).reshape(-1, 1)
    y = np.sin(X).ravel() + rng.normal(0, 0.05, size=n)
    return X, y


class TestMLPStructure(unittest.TestCase):
    def test_paper_dimension_count(self):
        # Section 4 worked example: (7,35,10,5,1) MLP has 701 parameters.
        mlp = MLPStructure([7, 35, 10, 5, 1])
        self.assertEqual(mlp.n_genes, 701)

    def test_flatten_unflatten_roundtrip(self):
        mlp = MLPStructure([3, 5, 2])
        space = mlp.build_parameter_space()
        genes = space.sample()
        weights, biases = mlp.unflatten(genes)
        self.assertEqual([w.shape for w in weights], [(3, 5), (5, 2)])
        self.assertEqual([b.shape for b in biases], [(5,), (2,)])
        np.testing.assert_allclose(mlp.flatten(weights, biases), genes)

    def test_forward_shape(self):
        mlp = MLPStructure([4, 6, 1])
        genes = mlp.build_parameter_space().sample()
        y = mlp.forward(genes, np.random.randn(10, 4))
        self.assertEqual(y.shape, (10, 1))


class TestBNNRegressionProblem(unittest.TestCase):
    def test_fitness_is_finite_and_nonnegative(self):
        X, y = make_regression_data()
        problem = BNNRegressionProblem([1, 8, 1], X, y)
        genes = problem.space.sample()
        loss = problem.fitness(genes)
        self.assertTrue(np.isfinite(loss))
        self.assertGreaterEqual(loss, 0.0)

    def test_perfect_fit_scores_near_zero(self):
        # A single linear layer can fit a linear function exactly.
        X = np.linspace(-2, 2, 20).reshape(-1, 1)
        y = 3.0 * X.ravel() + 1.0
        problem = BNNRegressionProblem([1, 1], X, y, output_activation=lambda z: z)
        # W=[[3]], b=[1] exactly reproduces y.
        genes = np.array([3.0, 1.0])
        self.assertLess(problem.fitness(genes), 1e-20)


class TestMutationOperators(unittest.TestCase):
    def setUp(self):
        np.random.seed(1)
        self.genomes = [np.full(3, float(i)) for i in range(8)]

    def test_rand_1_shape(self):
        mutant = mutate_rand_1(self.genomes, 0, F=0.5)
        self.assertEqual(mutant.shape, (3,))

    def test_rand_2_shape(self):
        mutant = mutate_rand_2(self.genomes, 0, F=0.5)
        self.assertEqual(mutant.shape, (3,))

    def test_best_1_uses_best_vector(self):
        best = np.full(3, 99.0)
        mutant = mutate_best_1(self.genomes, 0, F=0.0, best=best)
        np.testing.assert_allclose(mutant, best)

    def test_best_2_uses_best_vector(self):
        best = np.full(3, 99.0)
        mutant = mutate_best_2(self.genomes, 0, F=0.0, best=best)
        np.testing.assert_allclose(mutant, best)

    def test_all_operators_registered(self):
        self.assertEqual(set(MUTATION_OPERATORS), {"rand/1", "rand/2", "best/1", "best/2"})

    def test_mcmc_noise_changes_vector(self):
        np.random.seed(2)
        mutant = np.zeros(100)
        noisy = add_mcmc_noise(mutant, sigma2=1e-4)
        self.assertFalse(np.allclose(mutant, noisy))
        self.assertLess(np.std(noisy), 0.1)  # small noise, not swamping


class TestStagnationTracker(unittest.TestCase):
    def test_detects_stagnation_on_flat_sequence(self):
        tracker = StagnationTracker(window=5)
        stagnant = False
        for _ in range(10):
            stagnant = tracker.update(10.0)  # never improves
        self.assertTrue(stagnant)
        self.assertGreater(tracker.iteration, 0)

    def test_no_stagnation_while_improving(self):
        tracker = StagnationTracker(window=5)
        stagnant = False
        for val in [10.0, 8.0, 6.0, 4.0, 2.0, 0.0]:
            stagnant = tracker.update(val)
        self.assertFalse(stagnant)
        self.assertEqual(tracker.iteration, 0)

    def test_sample_hyperparameters_in_range(self):
        for _ in range(20):
            F, CR, operator = sample_hyperparameters()
            self.assertGreaterEqual(F, 0.10)
            self.assertLessEqual(F, 0.90)
            self.assertGreaterEqual(CR, 0.5)
            self.assertLessEqual(CR, 0.9)
            self.assertIn(operator, MUTATION_OPERATORS)


class TestRefinement(unittest.TestCase):
    def test_svd_and_local_search_never_worsen_population_best(self):
        np.random.seed(3)
        X, y = make_regression_data(n=20)
        problem = BNNRegressionProblem([1, 4, 1], X, y)
        solver = DEMCMCSolver(problem, options={"nPops": 8, "nGen": 1})
        solver.initialize()
        before = solver.population.population_sorted[0][1]

        svd_refine(solver.population, problem, before)
        local_search_refine(solver.population, problem,
                             solver.population.population_sorted[0][1])
        after = solver.population.population_sorted[0][1]
        self.assertLessEqual(after, before)


class TestDEMCMCSolver(unittest.TestCase):
    def test_requires_minimum_population(self):
        X, y = make_regression_data()
        problem = BNNRegressionProblem([1, 4, 1], X, y)
        with self.assertRaises(ValueError):
            DEMCMCSolver(problem, options={"nPops": 4, "nGen": 5})

    def test_rejects_unknown_mutation_operator(self):
        X, y = make_regression_data()
        problem = BNNRegressionProblem([1, 4, 1], X, y)
        with self.assertRaises(ValueError):
            DEMCMCSolver(problem, options={"nPops": 6, "nGen": 5, "mutation_operator": "bogus"})

    def test_converges_on_sine_regression(self):
        np.random.seed(4)
        X, y = make_regression_data(n=40)
        problem = BNNRegressionProblem([1, 8, 1], X, y)
        solver = DEMCMCSolver(problem, options={"nPops": 20, "nGen": 40})
        result = solver.run()
        self.assertLess(result.historyBest[-1], result.historyBest[0])

    def test_posterior_sample_count_single_chain(self):
        np.random.seed(5)
        X, y = make_regression_data(n=20)
        problem = BNNRegressionProblem([1, 4, 1], X, y)
        solver = DEMCMCSolver(problem, options={
            "nPops": 8, "nGen": 20, "burn_in": 12, "num_chains": 1,
        })
        solver.run()
        self.assertEqual(len(solver.posterior), 20 - 12)

    def test_posterior_sample_count_multi_chain(self):
        np.random.seed(6)
        X, y = make_regression_data(n=20)
        problem = BNNRegressionProblem([1, 4, 1], X, y)
        solver = DEMCMCSolver(problem, options={
            "nPops": 8, "nGen": 20, "burn_in": 12, "num_chains": 3,
        })
        solver.run()
        self.assertEqual(len(solver.posterior), (20 - 12) * 3)

    def test_hyper_mutation_and_refinement_do_not_crash(self):
        np.random.seed(7)
        X, y = make_regression_data(n=20)
        problem = BNNRegressionProblem([1, 4, 1], X, y)
        solver = DEMCMCSolver(problem, options={
            "nPops": 8, "nGen": 25, "burn_in": 10,
            "hyper_mutation": True, "stagnation_window": 3,
            "svd": True, "local_search": True, "refine_every": 5,
        })
        result = solver.run()
        self.assertEqual(len(result.historyBest), 25)

    def test_all_mutation_operators_run(self):
        X, y = make_regression_data(n=20)
        for operator in MUTATION_OPERATORS:
            with self.subTest(operator=operator):
                np.random.seed(8)
                problem = BNNRegressionProblem([1, 4, 1], X, y)
                solver = DEMCMCSolver(problem, options={
                    "nPops": 8, "nGen": 5, "mutation_operator": operator,
                })
                result = solver.run()
                self.assertEqual(len(result.historyBest), 5)


class TestPosteriorResult(unittest.TestCase):
    def test_predict_and_credible_interval_bracket_truth_reasonably(self):
        np.random.seed(9)
        X, y = make_regression_data(n=40)
        problem = BNNRegressionProblem([1, 8, 1], X, y)
        solver = DEMCMCSolver(problem, options={
            "nPops": 25, "nGen": 80, "burn_in": 30, "num_chains": 2,
        })
        solver.run()

        mean_pred, preds = solver.posterior.predict(problem)
        self.assertEqual(mean_pred.shape, (40, 1))
        lo, hi = solver.posterior.credible_interval(problem)
        self.assertTrue(np.all(lo <= hi))

        true_y = np.sin(X).ravel()
        mse = np.mean((mean_pred.ravel() - true_y) ** 2)
        self.assertLess(mse, 1.0)

    def test_mean_and_mode_genes_shapes(self):
        result = PosteriorResult()
        result.samples = [np.array([1.0, 2.0]), np.array([1.2, 1.8]), np.array([0.8, 2.2])]
        self.assertEqual(result.mean_genes().shape, (2,))
        self.assertEqual(result.mode_genes(bins=3).shape, (2,))

    def test_empty_posterior_raises(self):
        result = PosteriorResult()
        with self.assertRaises(ValueError):
            result.as_array()


class TestRegistry(unittest.TestCase):
    def test_lookup_de_mcmc(self):
        self.assertIs(get_solver("DE_MCMC"), DEMCMCSolver)
        self.assertIs(get_solver("de_mcmc"), DEMCMCSolver)


if __name__ == "__main__":
    unittest.main()
