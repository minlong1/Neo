import unittest

from Solvers import get_solver
from Solvers.core import OptimizationProblem

from PhysicsModules.XPS.xps_neo.problem import XPSProblem


class TestXPSScaffold(unittest.TestCase):
    def test_implements_problem_contract(self):
        problem = XPSProblem()
        self.assertIsInstance(problem, OptimizationProblem)
        self.assertGreaterEqual(problem.space.n_genes, 1)

    def test_pluggable_into_solver(self):
        problem = XPSProblem()
        solver = get_solver("GA")(problem, options={"nPops": 4, "nGen": 2})
        self.assertIsNotNone(solver.population)

    def test_fitness_not_implemented(self):
        with self.assertRaises(NotImplementedError):
            XPSProblem().fitness([0.0])


if __name__ == "__main__":
    unittest.main()
