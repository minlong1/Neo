"""
Tests for the AstroNeo OptimizationProblem plumbing.

Most of this module needs no PyXspec: build_parameter_space() and
genes_to_xspec_params() are pure numpy/dataclass code, independent of any
particular ModelSpec's real XSPEC parameter indices (those only exist once
a model is actually built). The end-to-end fitness tests do need PyXspec
(only installable via a HEASOFT build, see the module README) plus a real
spectrum, neither of which is available in CI -- guarded to skip unless
both PyXspec is importable and ASTRO_NEO_TEST_DATA_DIR points at a
directory with a usable spectrum file (ASTRO_NEO_TEST_DATA_FILE).
"""

import os
import unittest

import numpy as np

from PhysicsModules.AstroNeo.astro_neo.model import (
    DEFAULT_MODEL_SPEC,
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
        self.assertEqual([space.name_of(i) for i in range(space.n_genes)], list(PARAM_NAMES))

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


class TestDefaultModelSpec(unittest.TestCase):
    """DEFAULT_MODEL_SPEC's internal consistency -- pure dataclass checks,
    no XSPEC needed (real index resolution is covered by the live tests
    below, against DEFAULT_MODEL_SPEC's own known historical indices)."""

    def test_every_free_param_has_a_range(self):
        for name in DEFAULT_MODEL_SPEC.free_params:
            self.assertIn(name, DEFAULT_MODEL_SPEC.param_ranges)

    def test_every_free_and_linked_param_has_a_fixed_default(self):
        # fixed_params seeds every parameter (even ones later overwritten by
        # a free-param setPars() call or a link) -- same construction order
        # as the original hand-written model.
        linked_targets = {target for target, _ in DEFAULT_MODEL_SPEC.links}
        for name in DEFAULT_MODEL_SPEC.free_params:
            self.assertIn(name, DEFAULT_MODEL_SPEC.fixed_params)
        for name in linked_targets:
            self.assertIn(name, DEFAULT_MODEL_SPEC.fixed_params)

    def test_link_sources_are_not_dangling(self):
        # Every link's source should be a real parameter this spec sets
        # somewhere (fixed_params or free_params), not a typo'd name.
        known = set(DEFAULT_MODEL_SPEC.fixed_params) | set(DEFAULT_MODEL_SPEC.free_params)
        for target, source in DEFAULT_MODEL_SPEC.links:
            self.assertIn(source, known, f"link source {source!r} for {target!r} is unknown")


class TestGenesToXspecParams(unittest.TestCase):
    def test_maps_in_declared_order(self):
        indices = [10, 20, 30, 40, 50]
        genes = np.arange(len(indices), dtype=float)
        mapped = genes_to_xspec_params(genes, indices)
        self.assertEqual(list(mapped.keys()), indices)
        self.assertEqual(list(mapped.values()), list(genes))

    def test_values_are_python_floats(self):
        indices = list(range(1, 16))
        genes = np.array([1, 2, 3] * 5, dtype=np.float32)
        mapped = genes_to_xspec_params(genes, indices)
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

    def test_free_param_indices_match_historical_indices(self):
        # DEFAULT_MODEL_SPEC's free_params are dotted names now, resolved to
        # XSPEC indices at problem-build time -- pin them against the exact
        # indices the pre-ModelSpec version of this module hardcoded
        # (verified correct then; this is a regression check on the dotted
        # names being the *same* parameters, not new/typo'd ones).
        from PhysicsModules.AstroNeo.astro_neo.problem import AstroNeoProblem

        problem = AstroNeoProblem(
            data_dir=os.environ["ASTRO_NEO_TEST_DATA_DIR"],
            data_file=os.environ["ASTRO_NEO_TEST_DATA_FILE"],
            acx2_path=os.environ.get("ASTRO_NEO_ACX2_PATH"),
        )
        self.assertEqual(
            problem.free_param_indices,
            [2, 3, 4, 7, 9, 10, 11, 12, 13, 19, 22, 23, 38, 41, 60],
        )

    def test_custom_model_spec_fits_a_different_model(self):
        # Proves the generalization: a model with a different expression,
        # different free parameters, no links, and no ACX2 dependency still
        # fits cleanly through the same AstroNeoProblem/ModelSpec machinery.
        from PhysicsModules.AstroNeo.astro_neo.model import ModelSpec
        from PhysicsModules.AstroNeo.astro_neo.problem import AstroNeoProblem

        custom_spec = ModelSpec(
            expr="TBabs*powerlaw",
            free_params=("powerlaw.PhoIndex", "powerlaw.norm"),
            param_ranges={
                "powerlaw.PhoIndex": (0.5, 3.0),
                "powerlaw.norm": (1e-5, 1e-2),
            },
            fixed_params={"TBabs.nH": 0.0279},
        )
        problem = AstroNeoProblem(
            data_dir=os.environ["ASTRO_NEO_TEST_DATA_DIR"],
            data_file=os.environ["ASTRO_NEO_TEST_DATA_FILE"],
            model_spec=custom_spec,
        )
        self.assertEqual(problem.space.n_genes, 2)
        score_a = problem.fitness(problem.sample_genes())
        score_b = problem.fitness(problem.sample_genes())
        for score in (score_a, score_b):
            self.assertTrue(np.isfinite(score))
            self.assertNotEqual(score, 0.0)
        self.assertNotEqual(score_a, score_b)


if __name__ == "__main__":
    unittest.main()
