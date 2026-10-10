import os
import tempfile
import unittest

import numpy as np

from Solvers import get_solver
from Solvers.core import Individual, OptimizationProblem, Population, RunState

from PhysicsModules.Raman.raman_neo.problem import (
    RamanProblem,
    band_params,
    build_parameter_space,
    load_bands,
)
from PhysicsModules.Raman.raman_neo.bandrange import build_bandrange_pars
from PhysicsModules.Raman.tests.synthetic import (
    TRUE_BANDS,
    make_spectrum,
    true_genes,
    write_band_file,
    write_spectrum,
)


class ProblemTestBase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.dat = write_spectrum(os.path.join(cls.tmp.name, "spectrum.dat"))
        cls.bands = write_band_file(os.path.join(cls.tmp.name, "bands.txt"))

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def make_problem(self, **kwargs):
        defaults = dict(
            data_file=self.dat,
            band_file=self.bands,
            band_list=[1, 2, 3, 4],
            xmin=150.0,
            xmax=1750.0,
            deltax=1.0,
            normalize=False,  # keep the synthetic truth's absolute amplitudes
        )
        defaults.update(kwargs)
        return RamanProblem(**defaults)


class TestProblemContract(ProblemTestBase):
    def test_implements_problem_contract(self):
        problem = self.make_problem()
        self.assertIsInstance(problem, OptimizationProblem)
        self.assertEqual(problem.space.n_genes, 1 + 3 * 4)
        self.assertEqual(problem.space.name_of(0), "shift")
        self.assertEqual(problem.space.name_of(1), "band0_amp")
        self.assertEqual(problem.band_names, ["Band1", "Band2", "Band3", "Band4"])

    def test_default_band_list_uses_every_row(self):
        problem = self.make_problem(band_list=None)
        self.assertEqual(problem.nband, 6)

    def test_band_out_of_range_rejected(self):
        with self.assertRaises(ValueError):
            self.make_problem(band_list=[1, 9])

    def test_build_parameter_space(self):
        space = build_parameter_space(build_bandrange_pars(2), np.linspace(-1, 1, 21))
        self.assertEqual(space.n_genes, 7)
        self.assertEqual(space.limits(0), (-1.0, 1.0))

    def test_band_params_layout(self):
        genes = [2.5, 0.1, 5.0, -1.0, 0.2, 6.0, 1.0]
        params = band_params(genes, 2)
        np.testing.assert_allclose(params, [[0.1, 2.5, 5.0, -1.0], [0.2, 2.5, 6.0, 1.0]])

    def test_multi_compound_bands(self):
        a = write_band_file(os.path.join(self.tmp.name, "a.txt"), decoys=False)
        b = write_band_file(os.path.join(self.tmp.name, "b.txt"), decoys=False)
        bands = load_bands([a, b], [[1, 2], [3, 4]], num_compounds=2)
        self.assertEqual([x.name for x in bands], ["Comp0Band1", "Comp0Band2", "Comp1Band3", "Comp1Band4"])
        problem = self.make_problem(band_file=[a, b], band_list=[[1, 2], [3, 4]], num_compounds=2)
        self.assertLess(problem.fitness(true_genes(4)), 1e-12)
        with self.assertRaises(ValueError):
            load_bands([a], [[1, 2], [3, 4]], num_compounds=2)


class TestFitness(ProblemTestBase):
    def test_truth_scores_near_zero(self):
        problem = self.make_problem()
        self.assertLess(problem.fitness(true_genes(4)), 1e-12)

    def test_decoy_bands_with_zero_amp_are_harmless(self):
        problem = self.make_problem(band_list=[1, 2, 3, 4, 5, 6])
        self.assertLess(problem.fitness(true_genes(6)), 1e-12)

    def test_wrong_params_score_worse(self):
        problem = self.make_problem()
        truth = problem.fitness(true_genes(4))
        off = true_genes(4)
        off[0] = 4.0  # global shift
        off[2] = 20.0  # band0 width
        self.assertGreater(problem.fitness(off), truth + 1e-3)

    def test_weighting_changes_score(self):
        off = true_genes(4)
        off[1] = 0.1
        plain = self.make_problem().fitness(off)
        weighted = self.make_problem(weight=1.0).fitness(off)
        self.assertGreater(weighted, plain)

    def test_model_matches_data_on_window(self):
        problem = self.make_problem()
        np.testing.assert_allclose(problem.model_total(true_genes(4)), problem.y, atol=1e-9)
        curves = problem.band_curves(true_genes(4))
        np.testing.assert_allclose(np.sum(curves, axis=0), problem.y, atol=1e-9)

    def test_band_optimize_flags_decoys(self):
        problem = self.make_problem(band_list=[1, 2, 3, 4, 5, 6])
        result = problem.band_optimize(true_genes(6), 0.05)
        self.assertEqual(result["optimize_bands"], ["Band1", "Band2", "Band3", "Band4"])
        self.assertAlmostEqual(sum(result["optimize_band_lists"].values()), 1.0)


class TestShiftSweep(ProblemTestBase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.shifted = write_spectrum(
            os.path.join(cls.tmp.name, "shifted.dat"), make_spectrum(shift=3.0)
        )

    def test_optimize_shift_recovers_offset(self):
        np.random.seed(3)
        problem = self.make_problem(data_file=self.shifted)
        state = RunState(n_gen=10)
        population = Population(problem, state, n_pops=10)
        population.initialize_populations()

        start = true_genes(4, shift=0.0)
        state.globBestInd = Individual(problem.space, start)
        state.globBestVal = problem.fitness(start)

        shift = problem.optimize_shift(state, population)
        self.assertAlmostEqual(shift, 3.0, delta=0.06)
        self.assertAlmostEqual(state.globBestInd.genes[0], shift)
        self.assertLess(state.globBestVal, 1e-3)
        for individual in population.population:
            self.assertEqual(individual.genes[0], shift)

    def test_second_half_pins_sampled_shift(self):
        problem = self.make_problem()
        problem.best_shift = 1.5
        problem.second_half = True
        self.assertEqual(problem.sample_genes()[0], 1.5)
        problem.shift_sweep = False
        samples = {problem.sample_genes()[0] for _ in range(20)}
        self.assertGreater(len(samples), 1)

    def test_hooks_run_at_midpoint_and_end(self):
        np.random.seed(4)
        problem = self.make_problem(data_file=self.shifted)
        solver = get_solver("GA")(problem, options={"nPops": 20, "nGen": 4})
        solver.run()
        self.assertFalse(problem.second_half)  # reset by on_run_end
        self.assertAlmostEqual(solver.state.globBestInd.genes[0], problem.best_shift)


class TestConvergence(ProblemTestBase):
    def test_ga_reduces_misfit_and_finds_centers(self):
        np.random.seed(11)
        problem = self.make_problem()
        solver = get_solver("GA_RECHENBERG")(
            problem, options={"nPops": 80, "nGen": 30, "mutChance": 0.2}
        )
        result = solver.run()
        self.assertLess(result.historyBest[-1], result.historyBest[0])
        # the final shift sweep (on_run_end) runs after BaseSolver.run's last
        # collect: it updates the best genes and state, and can only improve
        score = problem.fitness(result.best_individual.genes)
        self.assertAlmostEqual(score, solver.state.globBestVal)
        self.assertLessEqual(score, result.best_value)

        params = problem.band_params(result.best_individual.genes)
        for band, (amp, shift, width, cshift), truth in zip(problem.bands, params, TRUE_BANDS):
            # the two strongest bands should land within a few cm^-1
            if truth[1] >= 0.8:
                self.assertAlmostEqual(band.center + shift + cshift, truth[0], delta=6.0)

    def test_every_solver_runs(self):
        for name in ("GA", "GA_RECHENBERG", "DE", "PSO", "CMA_ES"):
            with self.subTest(solver=name):
                np.random.seed(5)
                problem = self.make_problem()
                result = get_solver(name)(problem, options={"nPops": 12, "nGen": 4}).run()
                self.assertEqual(len(result.historyBest), 4)
                self.assertTrue(np.isfinite(result.best_value))


if __name__ == "__main__":
    unittest.main()
