"""
Tests for the AstroNeo OptimizationProblem plumbing.

Most of this module needs no PyXspec: build_parameter_space() and
genes_to_xspec_params() are pure numpy. The end-to-end fitness test does
need PyXspec (only installable via a HEASOFT build, see the module README)
plus a real spectrum, neither of which is available in CI -- it's guarded
to skip unless both PyXspec is importable and ASTRO_NEO_TEST_DATA_DIR points
at a directory with a usable spectrum file (ASTRO_NEO_TEST_DATA_FILE).
"""

import os
import unittest

import numpy as np

from PhysicsModules.AstroNeo.astro_neo.model import (
    FREE_PARAM_INDICES,
    PARAM_NAMES,
    PARAM_RANGES,
    genes_to_xspec_params,
)
from PhysicsModules.AstroNeo.astro_neo.problem import build_parameter_space


def _xspec_available():
    try:
        import xspec  # noqa: F401
    except ImportError:
        return False
    return True


class TestParameterSpace(unittest.TestCase):
    def test_gene_count_and_order(self):
        space = build_parameter_space()
        self.assertEqual(space.n_genes, len(PARAM_NAMES))
        self.assertEqual([space.name_of(i) for i in range(space.n_genes)], PARAM_NAMES)

    def test_gene_bounds_match_param_ranges(self):
        # GeneRange.from_bounds materializes np.arange(low, high, step), which
        # excludes the right endpoint -- gene_high sits one step short of
        # `high`, not exactly at it.
        space = build_parameter_space()
        for i, name in enumerate(PARAM_NAMES):
            low, high = PARAM_RANGES[name]
            step = (high - low) / 100_000
            gene_low, gene_high = space.limits(i)
            self.assertAlmostEqual(gene_low, low, places=6)
            self.assertLessEqual(abs(gene_high - high), step * 1.001)

    def test_sample_within_bounds(self):
        space = build_parameter_space()
        for _ in range(20):
            genes = space.sample()
            for i in range(space.n_genes):
                low, high = space.limits(i)
                self.assertGreaterEqual(genes[i], low)
                self.assertLessEqual(genes[i], high)


class TestGenesToXspecParams(unittest.TestCase):
    def test_maps_in_declared_order(self):
        genes = np.arange(len(PARAM_NAMES), dtype=float)
        mapped = genes_to_xspec_params(genes)
        self.assertEqual(list(mapped.keys()), FREE_PARAM_INDICES)
        self.assertEqual(list(mapped.values()), list(genes))

    def test_values_are_python_floats(self):
        genes = np.array([1, 2, 3] * 5, dtype=np.float32)
        mapped = genes_to_xspec_params(genes)
        self.assertTrue(all(isinstance(v, float) for v in mapped.values()))


@unittest.skipUnless(_xspec_available(), "PyXspec not importable (see README's manual install)")
@unittest.skipUnless(
    os.environ.get("ASTRO_NEO_TEST_DATA_DIR") and os.environ.get("ASTRO_NEO_TEST_DATA_FILE"),
    "set ASTRO_NEO_TEST_DATA_DIR/ASTRO_NEO_TEST_DATA_FILE to a real spectrum to run this",
)
class TestAstroNeoProblemLive(unittest.TestCase):
    """Only runs when explicitly pointed at real PyXspec + spectrum data;
    not part of the default/CI test run. See README "Testing" section."""

    def test_fitness_is_finite_nonzero_and_sensitive_to_genes(self):
        # A statistic stuck at exactly 0.0 (regardless of genes) means the
        # spectrum ended up with zero noticed channels -- e.g. the
        # ignore()/Plot.xAxis ordering bug this test caught once already
        # (Plot.xAxis must be set *before* .ignore(), since the ignore
        # string is parsed in whatever unit xAxis is at call time).
        from PhysicsModules.AstroNeo.astro_neo.problem import AstroNeoProblem

        problem = AstroNeoProblem(
            data_dir=os.environ["ASTRO_NEO_TEST_DATA_DIR"],
            data_file=os.environ["ASTRO_NEO_TEST_DATA_FILE"],
            acx2_path=os.environ.get("ASTRO_NEO_ACX2_PATH"),
        )
        score_a = problem.fitness(problem.sample_genes())
        score_b = problem.fitness(problem.sample_genes())
        for score in (score_a, score_b):
            self.assertTrue(np.isfinite(score))
            self.assertNotEqual(score, 0.0)
        self.assertNotEqual(score_a, score_b)


if __name__ == "__main__":
    unittest.main()
