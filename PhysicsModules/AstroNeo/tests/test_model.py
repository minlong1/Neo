"""
Tests for model.py's dotted-name parameter path resolution (_get_param/
_set_param) and the ModelSpec dataclass. The path resolvers are exercised
here with a plain fake object tree instead of a real xspec.Model -- they're
just a getattr/setattr chain walk, the mechanism the whole ModelSpec
generalization depends on, and worth pinning independently of whether
PyXspec is installed. Real xspec.Parameter behavior (`.index`, `.frozen`,
`.link`) is covered by test_problem.py's guarded live tests instead.
"""

import unittest

from PhysicsModules.AstroNeo.astro_neo.model import ModelSpec, _get_param, _set_param


class _Namespace:
    """Minimal stand-in for xspec's Component/Model attribute access --
    just plain nested objects, enough to exercise the dotted-path walk."""

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)


def _fake_model():
    return _Namespace(
        powerlaw=_Namespace(PhoIndex=1.5, norm=0.01),
        TBabs=_Namespace(nH=0.02),
        nested=_Namespace(sub=_Namespace(deep=42)),
    )


class TestGetParam(unittest.TestCase):
    def test_two_level_path(self):
        model = _fake_model()
        self.assertEqual(_get_param(model, "powerlaw.PhoIndex"), 1.5)
        self.assertEqual(_get_param(model, "powerlaw.norm"), 0.01)
        self.assertEqual(_get_param(model, "TBabs.nH"), 0.02)

    def test_deeper_path(self):
        model = _fake_model()
        self.assertEqual(_get_param(model, "nested.sub.deep"), 42)

    def test_unknown_component_raises_attribute_error(self):
        model = _fake_model()
        with self.assertRaises(AttributeError):
            _get_param(model, "doesnotexist.param")

    def test_unknown_parameter_raises_attribute_error(self):
        model = _fake_model()
        with self.assertRaises(AttributeError):
            _get_param(model, "powerlaw.doesnotexist")


class TestSetParam(unittest.TestCase):
    def test_two_level_path(self):
        model = _fake_model()
        _set_param(model, "powerlaw.PhoIndex", 2.7)
        self.assertEqual(model.powerlaw.PhoIndex, 2.7)

    def test_deeper_path(self):
        model = _fake_model()
        _set_param(model, "nested.sub.deep", 99)
        self.assertEqual(model.nested.sub.deep, 99)

    def test_does_not_disturb_sibling_attributes(self):
        model = _fake_model()
        _set_param(model, "powerlaw.PhoIndex", 9.9)
        self.assertEqual(model.powerlaw.norm, 0.01)
        self.assertEqual(model.TBabs.nH, 0.02)

    def test_get_after_set_round_trips(self):
        model = _fake_model()
        _set_param(model, "powerlaw.norm", 0.5)
        self.assertEqual(_get_param(model, "powerlaw.norm"), 0.5)


class TestModelSpecConstruction(unittest.TestCase):
    def _minimal_spec(self, **overrides):
        kwargs = dict(
            expr="powerlaw",
            free_params=("powerlaw.PhoIndex",),
            param_ranges={"powerlaw.PhoIndex": (0.0, 1.0)},
        )
        kwargs.update(overrides)
        return ModelSpec(**kwargs)

    def test_optional_fields_default_empty(self):
        spec = self._minimal_spec()
        self.assertEqual(spec.fixed_params, {})
        self.assertEqual(spec.unfrozen, ())
        self.assertEqual(spec.links, ())

    def test_fields_round_trip(self):
        spec = self._minimal_spec(
            fixed_params={"powerlaw.norm": 1.0},
            unfrozen=("powerlaw.norm",),
            links=(("a.b", "c.d"),),
        )
        self.assertEqual(spec.fixed_params, {"powerlaw.norm": 1.0})
        self.assertEqual(spec.unfrozen, ("powerlaw.norm",))
        self.assertEqual(spec.links, (("a.b", "c.d"),))

    def test_is_immutable(self):
        spec = self._minimal_spec()
        with self.assertRaises(Exception):
            spec.expr = "changed"


if __name__ == "__main__":
    unittest.main()
