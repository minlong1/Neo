"""
Tests for the AstroNeo runner's non-PyXspec surface: reading input (from a
dict or a real .ini file) and initializing output files. astro_setup()/
run()/plot() all construct an AstroNeoProblem and therefore need PyXspec;
those are covered by test_problem.py's guarded live tests instead.
"""

import os
import tempfile
import textwrap
import unittest
from pathlib import Path

from PhysicsModules.AstroNeo.astro_neo.astro_neo import AstroNeo

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


class TestAstroReadFromDict(unittest.TestCase):
    def test_sets_input_parameters(self):
        astro = AstroNeo()
        pars = {"data_dir": "x", "data_file": "y", "output_file": "z"}
        astro.astro_read(input_parameters=pars)
        self.assertEqual(astro.input_parameters, pars)

    def test_raises_without_filepath_or_dict(self):
        astro = AstroNeo()
        with self.assertRaises(ValueError):
            astro.astro_read()

    def test_explicit_dict_overrides_filepath_when_both_given(self):
        # astro_read applies filepath first, then input_parameters -- an
        # explicit dict always wins if both are passed.
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.ini"
            path.write_text(textwrap.dedent(MINIMAL_INI))
            astro = AstroNeo()
            override = {"data_dir": "override", "data_file": "y", "output_file": "z"}
            astro.astro_read(filepath=str(path), input_parameters=override)
        self.assertEqual(astro.input_parameters, override)


class TestAstroReadFromFile(unittest.TestCase):
    def test_reads_and_validates_ini_relative_to_cwd(self):
        with tempfile.TemporaryDirectory() as tmp:
            old_cwd = os.getcwd()
            os.chdir(tmp)
            try:
                Path("test.ini").write_text(textwrap.dedent(MINIMAL_INI))
                astro = AstroNeo()
                astro.astro_read(filepath="test.ini")
            finally:
                os.chdir(old_cwd)
        self.assertEqual(astro.input_parameters["data_file"], "spectrum.fits")
        self.assertEqual(astro.input_parameters["nPops"], 20)
        self.assertEqual(astro.input_parameters["solOpt"], 2)  # DE default

    def test_reads_ini_by_absolute_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.ini"
            path.write_text(textwrap.dedent(MINIMAL_INI))
            astro = AstroNeo()
            astro.astro_read(filepath=str(path))
        self.assertEqual(astro.input_parameters["data_dir"], "/some/dir")


class TestInitializeOutputs(unittest.TestCase):
    def test_creates_output_data_and_log_paths(self):
        with tempfile.TemporaryDirectory() as tmp:
            astro = AstroNeo()
            output_file = os.path.join(tmp, "run.csv")
            astro._initialize_outputs(output_file)

            self.assertTrue(os.path.exists(output_file))
            self.assertEqual(astro.data_path, os.path.join(tmp, "run_data.csv"))
            self.assertTrue(os.path.exists(os.path.join(tmp, "run.log")))

    def test_output_file_starts_with_header_row(self):
        with tempfile.TemporaryDirectory() as tmp:
            astro = AstroNeo()
            output_file = os.path.join(tmp, "run.csv")
            astro._initialize_outputs(output_file)
            content = Path(output_file).read_text()
        self.assertEqual(content, "Gen,TPS,CURRFIT,BESTFIT\n")

    def test_rerunning_clears_previous_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_file = os.path.join(tmp, "run.csv")
            AstroNeo()._initialize_outputs(output_file)
            with open(output_file, "a") as f:
                f.write("1,0.5,100.0,100.0\n")

            AstroNeo()._initialize_outputs(output_file)
            content = Path(output_file).read_text()
        # the stale data row from the "previous run" must be gone
        self.assertEqual(content, "Gen,TPS,CURRFIT,BESTFIT\n")

    def test_creates_nested_output_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            output_file = os.path.join(tmp, "nested", "dir", "run.csv")
            AstroNeo()._initialize_outputs(output_file)
            self.assertTrue(os.path.exists(output_file))


class TestCheckIfExists(unittest.TestCase):
    def test_removes_existing_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "f.txt")
            Path(path).write_text("stale")
            AstroNeo._check_if_exists(path)
            self.assertFalse(os.path.exists(path))

    def test_creates_parent_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "a", "b", "f.txt")
            AstroNeo._check_if_exists(path)
            self.assertTrue(os.path.isdir(os.path.join(tmp, "a", "b")))

    def test_no_error_when_file_absent(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "never_existed.txt")
            AstroNeo._check_if_exists(path)  # must not raise


if __name__ == "__main__":
    unittest.main()
