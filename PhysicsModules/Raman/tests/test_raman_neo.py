import csv
import os
import tempfile
import unittest

import numpy as np

from PhysicsModules.Raman.raman_neo.raman import RamanNeo
from PhysicsModules.Raman.tests.synthetic import (
    make_spectrum,
    write_band_file,
    write_spectrum,
)

INI_TEMPLATE = """\
[Inputs]
data_file = {data_file}
output_file = {output_file}
band_file = {band_file}

[Populations]
population = 30
num_gen = 6
best_sample = 9
lucky_few = 6

[Solver]
solver_options = {sol_opt}

[Mutations]
chance_of_mutation = 20
original_chance_of_mutation = 20
chance_of_mutation_e0 = 20
mutated_options = {mut_opt}

[Bands]
band_list = 1,2,3,4,5,6
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
"""


class TestRamanNeoEndToEnd(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        data = make_spectrum(baseline=lambda x: 0.05 + 5e-5 * x, noise=0.005, seed=1)
        self.dat = write_spectrum(os.path.join(self.tmp.name, "spectrum.dat"), data)
        self.bands = write_band_file(os.path.join(self.tmp.name, "bands.txt"))

    def tearDown(self):
        self.tmp.cleanup()

    def run_from_ini(self, sol_opt=1, mut_opt=1):
        output_file = os.path.join(self.tmp.name, "out", "run.csv")
        ini = os.path.join(self.tmp.name, "run.ini")
        with open(ini, "w") as f:
            f.write(
                INI_TEMPLATE.format(
                    data_file=self.dat,
                    band_file=self.bands,
                    output_file=output_file,
                    sol_opt=sol_opt,
                    mut_opt=mut_opt,
                )
            )

        raman = RamanNeo(verbose_lvl=0)
        raman.raman_read(filepath=ini)
        raman.raman_setup()
        result = raman.run()
        return raman, result, output_file

    def test_full_run_writes_outputs(self):
        np.random.seed(21)
        raman, result, output_file = self.run_from_ini()

        self.assertEqual(len(result.historyBest), 6)
        self.assertLessEqual(result.historyBest[-1], result.historyBest[0])
        self.assertAlmostEqual(
            raman.problem.fitness(result.best_individual.genes), result.best_value
        )

        with open(output_file) as f:
            lines = f.read().strip().splitlines()
        self.assertTrue(lines[0].startswith("Gen,TPS,"))
        self.assertEqual(len(lines), 1 + 6)  # header + one row per generation

        base = os.path.splitext(output_file)[0]
        with open(base + "_data.csv") as f:
            blocks = f.read().strip().split("#################################")
        self.assertEqual(len([b for b in blocks if b.strip()]), 6)
        self.assertEqual(len(blocks[0].strip().splitlines()), 6)  # one row per band
        self.assertTrue(os.path.exists(base + ".log"))

        with open(base + "_peaks.csv") as f:
            rows = list(csv.DictReader(f))
        self.assertEqual([r["band"] for r in rows], [f"Band{i}" for i in range(1, 7)])

        self.assertIn("optimize_bands", raman.band_optimize_result)
        self.assertIn("Raman Neo fit report", raman.analysis.report())

    def test_solver_and_mutation_options_run(self):
        for sol_opt, mut_opt in ((0, 0), (0, 1), (0, 3), (0, 4), (1, 1), (2, 1), ("PSO", 1)):
            with self.subTest(sol_opt=sol_opt, mut_opt=mut_opt):
                np.random.seed(30)
                raman, result, _ = self.run_from_ini(sol_opt=sol_opt, mut_opt=mut_opt)
                self.assertEqual(len(result.historyBest), 6)
                self.assertTrue(np.isfinite(result.best_value))

    def test_dict_inputs(self):
        np.random.seed(7)
        raman = RamanNeo(verbose_lvl=0)
        raman.raman_read(
            input_parameters={
                "data_file": self.dat,
                "band_file": self.bands,
                "band_list": [1, 2, 3, 4],
                "output_file": os.path.join(self.tmp.name, "dict", "run.csv"),
                "nPops": 20,
                "nGen": 3,
                "xmin": 150,
                "xmax": 1750,
                "baseline_order": 1,
            }
        )
        raman.raman_setup()
        result = raman.run()
        self.assertEqual(raman.problem.nband, 4)
        self.assertEqual(len(result.historyBest), 3)
        self.assertIsNone(raman.band_optimize_result)

    def test_missing_inputs_raise(self):
        with self.assertRaises(ValueError):
            RamanNeo(verbose_lvl=0).raman_read()


if __name__ == "__main__":
    unittest.main()
