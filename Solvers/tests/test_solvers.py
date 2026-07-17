import unittest

import numpy as np

from Solvers import DESolver, GARechenbergSolver, GASolver, PSOSolver, get_solver
from Solvers.tests.toy_problem import QuadraticProblem


class TestRegistry(unittest.TestCase):
    def test_lookup_by_name(self):
        self.assertIs(get_solver("GA"), GASolver)
        self.assertIs(get_solver("ga"), GASolver)
        self.assertIs(get_solver("de"), DESolver)
        self.assertIs(get_solver("pso"), PSOSolver)

    def test_lookup_by_numeric_id(self):
        self.assertIs(get_solver(0), GASolver)
        self.assertIs(get_solver(1), GARechenbergSolver)
        self.assertIs(get_solver(2), DESolver)

    def test_unknown_solver(self):
        with self.assertRaises(ValueError):
            get_solver("simulated_quantum_annealing")


class TestGASolverConvergence(unittest.TestCase):
    def test_converges_on_quadratic(self):
        np.random.seed(42)
        problem = QuadraticProblem()
        solver = GASolver(
            problem,
            options={"nPops": 60, "nGen": 30, "mutChance": 0.3},
        )
        result = solver.run()

        initial = result.historyBest[0]
        final = result.historyBest[-1]
        self.assertLess(final, initial)
        # Grid resolution is 0.01 over 4 genes; GA should get close to 0
        self.assertLess(final, 0.05)
        self.assertEqual(len(result.historyBest), 30)

    def test_history_monotonic_nonincreasing(self):
        np.random.seed(1)
        problem = QuadraticProblem()
        result = GASolver(problem, options={"nPops": 30, "nGen": 15}).run()
        for prev, curr in zip(result.historyBest, result.historyBest[1:]):
            self.assertLessEqual(curr, prev)

    def test_hooks_fire(self):
        np.random.seed(2)
        problem = QuadraticProblem()
        GASolver(problem, options={"nPops": 10, "nGen": 5}).run()
        self.assertEqual(problem.generation_end_calls, 5)
        self.assertEqual(problem.run_end_calls, 1)


class TestGARechenberg(unittest.TestCase):
    def test_runs_and_adapts(self):
        np.random.seed(3)
        problem = QuadraticProblem()
        solver = GARechenbergSolver(
            problem, options={"nPops": 30, "nGen": 25, "mutChance": 0.3}
        )
        result = solver.run()
        self.assertEqual(len(result.historyBest), 25)
        self.assertGreaterEqual(solver.mutator.mut_chance, 0.0)
        self.assertLessEqual(solver.mutator.mut_chance, 1.0)


class TestDESolver(unittest.TestCase):
    def test_converges_on_quadratic(self):
        np.random.seed(4)
        problem = QuadraticProblem()
        solver = DESolver(problem, options={"nPops": 40, "nGen": 40})
        result = solver.run()

        initial = result.historyBest[0]
        final = result.historyBest[-1]
        self.assertLess(final, initial)
        self.assertLess(final, 0.05)
        self.assertEqual(len(result.historyBest), 40)

    def test_history_monotonic_nonincreasing(self):
        np.random.seed(5)
        problem = QuadraticProblem()
        result = DESolver(problem, options={"nPops": 20, "nGen": 15}).run()
        for prev, curr in zip(result.historyBest, result.historyBest[1:]):
            self.assertLessEqual(curr, prev)

    def test_hooks_fire(self):
        np.random.seed(6)
        problem = QuadraticProblem()
        DESolver(problem, options={"nPops": 10, "nGen": 5}).run()
        self.assertEqual(problem.generation_end_calls, 5)
        self.assertEqual(problem.run_end_calls, 1)

    def test_trial_genes_stay_in_bounds(self):
        np.random.seed(7)
        problem = QuadraticProblem()
        solver = DESolver(problem, options={"nPops": 10, "nGen": 10, "F": 1.5})
        result = solver.run()
        for i in range(len(result.best_individual.genes)):
            low, high = problem.space.limits(i)
            self.assertGreaterEqual(result.best_individual.genes[i], low)
            self.assertLessEqual(result.best_individual.genes[i], high)

    def test_requires_at_least_four_individuals(self):
        problem = QuadraticProblem()
        with self.assertRaises(ValueError):
            DESolver(problem, options={"nPops": 3, "nGen": 5})


class TestPSOSolver(unittest.TestCase):
    def test_converges_on_quadratic(self):
        np.random.seed(8)
        problem = QuadraticProblem()
        solver = PSOSolver(problem, options={"nPops": 40, "nGen": 40})
        result = solver.run()

        initial = result.historyBest[0]
        final = result.historyBest[-1]
        self.assertLess(final, initial)
        self.assertLess(final, 0.05)
        self.assertEqual(len(result.historyBest), 40)

    def test_history_monotonic_nonincreasing(self):
        np.random.seed(9)
        problem = QuadraticProblem()
        result = PSOSolver(problem, options={"nPops": 20, "nGen": 15}).run()
        for prev, curr in zip(result.historyBest, result.historyBest[1:]):
            self.assertLessEqual(curr, prev)

    def test_particle_genes_stay_in_bounds(self):
        np.random.seed(10)
        problem = QuadraticProblem()
        solver = PSOSolver(problem, options={"nPops": 10, "nGen": 10})
        result = solver.run()
        for i in range(len(result.best_individual.genes)):
            low, high = problem.space.limits(i)
            self.assertGreaterEqual(result.best_individual.genes[i], low)
            self.assertLessEqual(result.best_individual.genes[i], high)

    def test_requires_at_least_two_individuals(self):
        problem = QuadraticProblem()
        with self.assertRaises(ValueError):
            PSOSolver(problem, options={"nPops": 1, "nGen": 5})


if __name__ == "__main__":
    unittest.main()
