import os
import tempfile
import unittest

import numpy as np

from PhysicsModules.NanoIndentation.nanoindentation_neo.nano_indent import NanoNeo
from PhysicsModules.NanoIndentation.tests.synthetic import write_csv

INI_TEMPLATE = """\
[Inputs]
data_file = {data_file}
output_file = {output_file}
data_cutoff = 0.1, 0.9

[Populations]
population = 40
num_gen = 5
best_sample = 30
lucky_few = 20

[Mutations]
chance_of_mutation = 30
original_chance_of_mutation = 30
mutated_options = {mut_opt}

[Paths]
npaths = 1
fits = OliverPharr
a_range = 1e-4, 3e-4, 1e-5
hf_range = 250.0, 350.0, 1.0
m_range = 1.3, 1.7, 0.01

[Outputs]
print_graph = False
num_output_paths = False
"""


class TestNanoNeoEndToEnd(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.csv = write_csv(os.path.join(self.tmp.name, "curve.csv"))

    def tearDown(self):
        self.tmp.cleanup()

    def run_from_ini(self, mut_opt=0):
        output_file = os.path.join(self.tmp.name, "out", "run.csv")
        ini = os.path.join(self.tmp.name, "run.ini")
        with open(ini, "w") as f:
            f.write(
                INI_TEMPLATE.format(
                    data_file=self.csv, output_file=output_file, mut_opt=mut_opt
                )
            )

        nano = NanoNeo()
        nano.nano_read(filepath=ini)
        nano.nano_setup()
        result = nano.run()
        return nano, result, output_file

    def test_full_run_writes_outputs(self):
        np.random.seed(21)
        nano, result, output_file = self.run_from_ini(mut_opt=0)

        self.assertEqual(len(result.historyBest), 5)
        self.assertLessEqual(result.historyBest[-1], result.historyBest[0])

        with open(output_file) as f:
            lines = f.read().strip().splitlines()
        self.assertTrue(lines[0].startswith("Gen,TPS,"))
        self.assertEqual(len(lines), 1 + 5)  # header + one row per generation

        data_file = os.path.splitext(output_file)[0] + "_data.csv"
        log_file = os.path.splitext(output_file)[0] + ".log"
        self.assertTrue(os.path.exists(data_file))
        self.assertTrue(os.path.exists(log_file))

    def test_all_mutation_options_run(self):
        for mut_opt in (0, 1, 2, 3):
            with self.subTest(mut_opt=mut_opt):
                np.random.seed(30 + mut_opt)
                _, result, _ = self.run_from_ini(mut_opt=mut_opt)
                self.assertEqual(len(result.historyBest), 5)


if __name__ == "__main__":
    unittest.main()
