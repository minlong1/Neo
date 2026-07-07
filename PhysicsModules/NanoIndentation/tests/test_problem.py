import os
import tempfile
import unittest

import numpy as np

from Solvers import get_solver
from Solvers.core import OptimizationProblem

from PhysicsModules.NanoIndentation.nanoindentation_neo.problem import (
    NanoIndentationProblem,
    build_parameter_space,
)
from PhysicsModules.NanoIndentation.tests.synthetic import (
    TRUE_A,
    TRUE_HF,
    TRUE_M,
    write_csv,
)


class ProblemTestBase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.csv = write_csv(os.path.join(cls.tmp.name, "curve.csv"))

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def make_problem(self, **kwargs):
        defaults = dict(
            data_file=self.csv,
            data_cutoff=(0.1, 0.9),
            npaths=1,
            # narrow grids bracketing the ground truth
            A_range=(1e-4, 3e-4, 1e-5),
            hf_range=(250.0, 350.0, 1.0),
            m_range=(1.3, 1.7, 0.01),
        )
        defaults.update(kwargs)
        return NanoIndentationProblem(**defaults)


class TestProblemContract(ProblemTestBase):
    def test_implements_problem_contract(self):
        problem = self.make_problem()
        self.assertIsInstance(problem, OptimizationProblem)
        self.assertEqual(problem.space.n_genes, 3)

    def test_multi_path_space(self):
        problem = self.make_problem(npaths=2)
        self.assertEqual(problem.space.n_genes, 6)

    def test_unknown_model_rejected(self):
        with self.assertRaises(ValueError):
            self.make_problem(fits="NotAModel")

    def test_build_parameter_space_defaults(self):
        space = build_parameter_space(1)
        self.assertEqual(space.n_genes, 3)


class TestFitness(ProblemTestBase):
    def test_truth_scores_near_zero(self):
        problem = self.make_problem()
        loss = problem.fitness([TRUE_A, TRUE_HF, TRUE_M])
        self.assertLess(loss, 1e-12)

    def test_wrong_params_score_worse(self):
        problem = self.make_problem()
        truth = problem.fitness([TRUE_A, TRUE_HF, TRUE_M])
        off = problem.fitness([TRUE_A * 1.5, TRUE_HF - 40.0, TRUE_M + 0.1])
        self.assertGreater(off, truth)

    def test_nan_region_scores_inf(self):
        # hf above the whole x slice makes (h - hf)^m NaN -> inf, not NaN
        problem = self.make_problem(hf_range=(250.0, 2000.0, 10.0))
        loss = problem.fitness([TRUE_A, 1900.0, 1.5])
        self.assertEqual(loss, np.inf)
        self.assertFalse(np.isnan(loss))


class TestConvergence(ProblemTestBase):
    def test_ga_converges_to_ground_truth(self):
        np.random.seed(11)
        problem = self.make_problem()
        solver = get_solver("GA")(
            problem, options={"nPops": 80, "nGen": 40, "mutChance": 0.3}
        )
        result = solver.run()

        self.assertLess(result.historyBest[-1], result.historyBest[0])
        A, hf, m = result.best_individual.genes
        self.assertAlmostEqual(hf, TRUE_HF, delta=10.0)
        self.assertAlmostEqual(m, TRUE_M, delta=0.1)

    def test_rechenberg_solver_runs(self):
        np.random.seed(12)
        problem = self.make_problem()
        result = get_solver("GA_RECHENBERG")(
            problem, options={"nPops": 40, "nGen": 10, "mutChance": 0.3}
        ).run()
        self.assertEqual(len(result.historyBest), 10)


if __name__ == "__main__":
    unittest.main()
