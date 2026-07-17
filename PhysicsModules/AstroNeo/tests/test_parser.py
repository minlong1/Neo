"""
Tests for Astro Neo's .ini parsing (parser.py). Pure Python, no PyXspec
needed -- these only exercise configparser plus the CheckKey/optional_var/
validate plumbing, not AstroNeoProblem.
"""

import tempfile
import textwrap
import unittest
from pathlib import Path

from PhysicsModules.AstroNeo.astro_neo.parser import (
    CheckKey,
    optional_var,
    parse_input_file,
    read_input_file,
)

MINIMAL_INI = """
[Inputs]
data_dir = /some/dir
data_file = spectrum.fits
output_file = out.csv

[Populations]
population = 20
num_gen = 20

[Mutations]
mutated_options = 4

[Outputs]
print_graph = False
"""

FULL_INI = """
[Inputs]
data_dir = /some/dir
data_file = spectrum.fits
bg_file = bg.fits
rsp_file = rsp.fits
acx2_path = /some/acx2
output_file = out.csv
xmin = 6.0
xmax = 25.0

[Populations]
population = 50
num_gen = 100
best_sample = 0.4
lucky_few = 0.1

[Mutations]
solver_type = 0
mutated_options = 1
crossover_options = 3
chance_of_mutation = 0.2
mutf = 0.7
mutcr = 0.8

[Outputs]
print_graph = True
distributed = 2
"""


def _write_ini(tmpdir, contents):
    path = Path(tmpdir) / "test.ini"
    path.write_text(textwrap.dedent(contents))
    return str(path)


class TestOptionalVar(unittest.TestCase):
    def test_returns_value_when_present(self):
        self.assertEqual(optional_var({"foo": "3.5"}, "foo", 1.0, float), 3.5)

    def test_returns_default_when_missing(self):
        self.assertEqual(optional_var({}, "foo", 1.0, float), 1.0)

    def test_int_conversion(self):
        self.assertEqual(optional_var({"n": "7"}, "n", 0, int), 7)

    def test_bool_true_and_false(self):
        self.assertTrue(optional_var({"flag": "True"}, "flag", False, bool))
        self.assertFalse(optional_var({"flag": "False"}, "flag", True, bool))

    def test_bool_missing_uses_default(self):
        self.assertFalse(optional_var({}, "flag", False, bool))
        self.assertTrue(optional_var({}, "flag", True, bool))

    def test_none_type_returns_raw_string_or_none(self):
        self.assertEqual(optional_var({"x": "raw"}, "x", None, None), "raw")
        self.assertIsNone(optional_var({}, "x", None, None))


class TestCheckKey(unittest.TestCase):
    def test_raises_on_missing_key(self):
        with self.assertRaises(KeyError):
            CheckKey({"a": 1}, ["a", "b"])

    def test_passes_when_all_present(self):
        CheckKey({"a": 1, "b": 2}, ["a", "b"])  # must not raise


class TestReadInputFile(unittest.TestCase):
    def test_missing_file_raises_file_not_found(self):
        with self.assertRaises(FileNotFoundError):
            read_input_file("/nonexistent/path/to/astro_neo_test.ini")

    def test_missing_section_raises_key_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = _write_ini(
                tmp,
                """
                [Inputs]
                data_dir = /x
                data_file = y.fits
                output_file = z.csv
                """,
            )
            with self.assertRaises(KeyError):
                read_input_file(path)

    def test_missing_required_key_in_present_section_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            # Inputs section present but missing data_file.
            path = _write_ini(
                tmp,
                """
                [Inputs]
                data_dir = /x
                output_file = z.csv

                [Populations]
                population = 20
                num_gen = 20

                [Mutations]
                mutated_options = 4

                [Outputs]
                print_graph = False
                """,
            )
            with self.assertRaises(KeyError):
                read_input_file(path)

    def test_reads_minimal_valid_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = _write_ini(tmp, MINIMAL_INI)
            config = read_input_file(path)
        self.assertEqual(config["Inputs"]["data_file"], "spectrum.fits")
        self.assertEqual(config["Populations"]["population"], "20")


class TestValidateInputFile(unittest.TestCase):
    def test_minimal_file_gets_defaults(self):
        with tempfile.TemporaryDirectory() as tmp:
            pars = parse_input_file(_write_ini(tmp, MINIMAL_INI))

        self.assertEqual(pars["data_dir"], "/some/dir")
        self.assertEqual(pars["data_file"], "spectrum.fits")
        self.assertIsNone(pars["bg_file"])
        self.assertIsNone(pars["rsp_file"])
        self.assertIsNone(pars["acx2_path"])
        self.assertEqual(pars["xmin"], 7.0)
        self.assertEqual(pars["xmax"], 30.0)
        self.assertEqual(pars["nPops"], 20)
        self.assertEqual(pars["nGen"], 20)
        self.assertEqual(pars["nBestSample"], 0.3)
        self.assertEqual(pars["nLuckSample"], 0.2)
        self.assertEqual(pars["solOpt"], 2)
        self.assertEqual(pars["mutOpt"], 4)
        self.assertEqual(pars["croOpt"], 6)
        self.assertEqual(pars["mutChance"], 0.3)
        self.assertEqual(pars["F"], 0.5)
        self.assertEqual(pars["CR"], 0.9)
        self.assertFalse(pars["printGraph"])
        self.assertEqual(pars["distributed"], 1)

    def test_full_file_overrides_every_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            pars = parse_input_file(_write_ini(tmp, FULL_INI))

        self.assertEqual(pars["bg_file"], "bg.fits")
        self.assertEqual(pars["rsp_file"], "rsp.fits")
        self.assertEqual(pars["acx2_path"], "/some/acx2")
        self.assertEqual(pars["xmin"], 6.0)
        self.assertEqual(pars["xmax"], 25.0)
        self.assertEqual(pars["nPops"], 50)
        self.assertEqual(pars["nGen"], 100)
        self.assertEqual(pars["nBestSample"], 0.4)
        self.assertEqual(pars["nLuckSample"], 0.1)
        self.assertEqual(pars["solOpt"], 0)
        self.assertEqual(pars["mutOpt"], 1)
        self.assertEqual(pars["croOpt"], 3)
        self.assertAlmostEqual(pars["mutChance"], 0.2)
        self.assertAlmostEqual(pars["F"], 0.7)
        self.assertAlmostEqual(pars["CR"], 0.8)
        self.assertTrue(pars["printGraph"])
        self.assertEqual(pars["distributed"], 2)

    def test_returned_types_are_correctly_coerced(self):
        with tempfile.TemporaryDirectory() as tmp:
            pars = parse_input_file(_write_ini(tmp, FULL_INI))

        self.assertIsInstance(pars["nPops"], int)
        self.assertIsInstance(pars["nGen"], int)
        self.assertIsInstance(pars["xmin"], float)
        self.assertIsInstance(pars["mutChance"], float)
        self.assertIsInstance(pars["printGraph"], bool)


if __name__ == "__main__":
    unittest.main()
