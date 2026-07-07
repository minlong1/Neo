import os
import tempfile
import unittest

from PhysicsModules.NanoIndentation.nanoindentation_neo.parser import (
    parse_input_file,
    read_input_file,
)

SAMPLE_INI = """\
[Inputs]
data_file = data/sample.csv
output_file = result/out.csv
data_cutoff = 0.2, 0.95

[Populations]
population = 1800
num_gen = 70
best_sample = 20
lucky_few = 10

[Mutations]
chance_of_mutation = 30
original_chance_of_mutation = 20
mutated_options = 0

[Paths]
npaths = 1
fits = OliverPharr
a_range = 0.01, 100.0, 0.01
hf_range = 270.0, 570.0, 0.1
m_range = 1.1, 1.7, 0.01

[Outputs]
print_graph = False
num_output_paths = False
"""


class TestParser(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.ini = os.path.join(self.tmp.name, "test.ini")
        with open(self.ini, "w") as f:
            f.write(SAMPLE_INI)

    def tearDown(self):
        self.tmp.cleanup()

    def test_read_sections(self):
        sections = read_input_file(self.ini)
        self.assertIn("Inputs", sections)
        self.assertIn("Populations", sections)

    def test_parse_values(self):
        pars = parse_input_file(self.ini)
        self.assertEqual(pars["data_file"], "data/sample.csv")
        self.assertEqual(pars["data_cutoff"], [0.2, 0.95])
        self.assertEqual(pars["nPops"], 1800)
        self.assertEqual(pars["nGen"], 70)
        self.assertAlmostEqual(pars["nBestSample"], 0.20)
        self.assertAlmostEqual(pars["nLuckSample"], 0.10)
        self.assertAlmostEqual(pars["mutChance"], 0.30)
        self.assertEqual(pars["mutOpt"], 0)
        self.assertEqual(pars["npaths"], 1)
        self.assertEqual(pars["fits"], "OliverPharr")
        self.assertEqual(pars["A_range"], [0.01, 100.0, 0.01])
        self.assertEqual(pars["hf_range"], [270.0, 570.0, 0.1])
        self.assertEqual(pars["m_range"], [1.1, 1.7, 0.01])
        self.assertFalse(pars["printGraph"])

    def test_missing_section_raises(self):
        bad = os.path.join(self.tmp.name, "bad.ini")
        with open(bad, "w") as f:
            f.write("[Inputs]\ndata_file = x\noutput_file = y\ndata_cutoff = 0.1,0.9\n")
        with self.assertRaises(KeyError):
            read_input_file(bad)

    def test_missing_file_raises(self):
        with self.assertRaises(FileNotFoundError):
            read_input_file(os.path.join(self.tmp.name, "nope.ini"))


if __name__ == "__main__":
    unittest.main()
