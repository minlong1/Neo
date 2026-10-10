import os
import tempfile
import unittest

from PhysicsModules.Raman.raman_neo.parser import (
    parse_input_file,
    read_input_file,
    solver_option,
    split_band_arr,
)

SAMPLE_INI = """\
[Inputs]
data_file = data/sample_raman.dat
output_file = result/out.csv
band_file = data/sample_bands.txt
num_compounds = 1

[Populations]
population = 100
num_gen = 40
best_sample = 30
lucky_few = 20

[Solver]
solver_options = 1

[Mutations]
chance_of_mutation = 20
original_chance_of_mutation = 20
chance_of_mutation_e0 = 20
mutated_options = 1
selection_options = 0
crossover_options = 0

[Bands]
band_list = 1,2,3,4,5,6
individual_band = False
band_optimize = True
band_optimize_percent = 0.05

[Spectrum]
xmin = 150
xmax = 1750
weight = 0.0
deltax = 1.0
baseline_order = 1

[Outputs]
print_graph = False
num_output_bands = True
steady_state_exit = False
"""


class TestParser(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.ini = self.write(SAMPLE_INI)

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, text, name="test.ini"):
        path = os.path.join(self.tmp.name, name)
        with open(path, "w") as f:
            f.write(text)
        return path

    def test_read_sections(self):
        sections = read_input_file(self.ini)
        for name in ("Inputs", "Populations", "Solver", "Mutations", "Bands", "Spectrum", "Outputs"):
            self.assertIn(name, sections)

    def test_parse_values(self):
        pars = parse_input_file(self.ini)
        self.assertEqual(pars["data_file"], "data/sample_raman.dat")
        self.assertEqual(pars["band_file"], "data/sample_bands.txt")
        self.assertIsNone(pars["bandrange_file"])
        self.assertEqual(pars["num_compounds"], 1)
        self.assertEqual(pars["nPops"], 100)
        self.assertEqual(pars["nGen"], 40)
        # best_sample / lucky_few are counts relative to the population
        self.assertAlmostEqual(pars["nBestSample"], 0.30)
        self.assertAlmostEqual(pars["nLuckSample"], 0.20)
        self.assertEqual(pars["solOpt"], 1)
        self.assertAlmostEqual(pars["mutChance"], 0.20)
        self.assertEqual(pars["mutOpt"], 1)
        self.assertEqual(pars["selOpt"], 0)
        self.assertEqual(pars["croOpt"], 0)
        self.assertEqual(pars["band_list"], [1, 2, 3, 4, 5, 6])
        self.assertTrue(pars["bandOptimize"])
        self.assertAlmostEqual(pars["bandOptimizePercent"], 0.05)
        self.assertEqual(pars["xmin"], 150.0)
        self.assertEqual(pars["xmax"], 1750.0)
        self.assertEqual(pars["baseline_order"], 1)
        self.assertTrue(pars["normalize"])
        self.assertEqual(pars["shift_range"], [])
        self.assertFalse(pars["printGraph"])
        self.assertFalse(pars["steadyState"])

    def test_optional_keys_default(self):
        text = SAMPLE_INI.replace("solver_options = 1\n", "").replace("baseline_order = 1\n", "")
        pars = parse_input_file(self.write(text, "min.ini"))
        self.assertEqual(pars["solOpt"], 0)
        self.assertEqual(pars["baseline_order"], 0)

    def test_shift_range_and_named_solver(self):
        text = SAMPLE_INI.replace("solver_options = 1", "solver_options = PSO").replace(
            "baseline_order = 1", "baseline_order = 1\nshift_range = -5, 5, 0.5"
        )
        pars = parse_input_file(self.write(text, "named.ini"))
        self.assertEqual(pars["solOpt"], "PSO")
        self.assertEqual(pars["shift_range"], [-5.0, 5.0, 0.5])

    def test_multi_compound(self):
        text = SAMPLE_INI.replace("num_compounds = 1", "num_compounds = 2").replace(
            "band_file = data/sample_bands.txt", "band_file = a.txt, b.txt"
        ).replace("band_list = 1,2,3,4,5,6", "band_list = [1,2][3,4]")
        pars = parse_input_file(self.write(text, "multi.ini"))
        self.assertEqual(pars["band_file"], ["a.txt", "b.txt"])
        self.assertEqual(pars["band_list"], [[1, 2], [3, 4]])

    def test_split_band_arr(self):
        self.assertEqual(split_band_arr("1, 2,3"), [1, 2, 3])
        self.assertEqual(split_band_arr("[1,2]"), [1, 2])
        self.assertEqual(split_band_arr("[1,2][3]", 2), [[1, 2], [3]])
        with self.assertRaises(ValueError):
            split_band_arr("[1,2][3]", 3)
        with self.assertRaises(ValueError):
            split_band_arr("[1,2", 1)

    def test_solver_option(self):
        self.assertEqual(solver_option("2"), 2)
        self.assertEqual(solver_option(" cma_es "), "cma_es")

    def test_missing_section_raises(self):
        bad = self.write("[Inputs]\ndata_file = x\noutput_file = y\nband_file = z\n", "bad.ini")
        with self.assertRaises(KeyError):
            read_input_file(bad)

    def test_missing_key_raises(self):
        bad = self.write(SAMPLE_INI.replace("band_list = 1,2,3,4,5,6\n", ""), "nokey.ini")
        with self.assertRaises(KeyError):
            read_input_file(bad)

    def test_missing_file_raises(self):
        with self.assertRaises(FileNotFoundError):
            read_input_file(os.path.join(self.tmp.name, "nope.ini"))


if __name__ == "__main__":
    unittest.main()
