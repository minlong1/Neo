import unittest

from PhysicsModules.EXAFS.exafs_neo.utils_mapping import (
    neocrossover_int2str, neocrossover_str2int,
    neomutator_int2str, neomutator_str2int,
    neoselector_int2str, neoselector_str2int,
    neosolver_int2str, neosolver_str2int,
)

CROSSOVER = [
    (0, "Uniform Crossover"),
    (1, "Single Point Crossover"),
    (2, "Dual Point Crossover"),
    (3, "Arithmetic Crossover"),
    (4, "Or Crossover"),
    (5, "Average Crossover"),
]

MUTATOR = [
    (0, "Mutate Per Individual"),
    (1, "Mutate Per Path"),
    (2, "Mutate Per Trait"),
    (3, "Mutate Metropolis"),
    (4, "Mutate Bounded Per Range"),
    (5, "Mutate Differential Evolution"),
]

SELECTOR = [
    (0, "Roulette Wheel"),
    (1, "Tournament"),
]

SOLVER = [
    (0, "Genetic Algorithm"),
    (1, "Genetic Algorithm with Rechenberg"),
    (2, "Differential Evolution"),
]


class TestUtilsMapping(unittest.TestCase):
    def test_neocrossover_mappings(self):
        for i, s in CROSSOVER:
            self.assertEqual(neocrossover_int2str(i), s)
            self.assertEqual(neocrossover_str2int(s), i)

    def test_neocrossover_roundtrip_and_unknowns(self):
        for i, _ in CROSSOVER:
            self.assertEqual(neocrossover_str2int(neocrossover_int2str(i)), i)
        self.assertEqual(neocrossover_int2str(999), "Unknown Crossover")
        self.assertEqual(neocrossover_int2str(None), "Unknown Crossover")
        self.assertEqual(neocrossover_str2int("nope"), -1)
        self.assertEqual(neocrossover_str2int(None), -1)

    def test_neomutator_mappings(self):
        for i, s in MUTATOR:
            self.assertEqual(neomutator_int2str(i), s)
            self.assertEqual(neomutator_str2int(s), i)

    def test_neomutator_roundtrip_and_unknowns(self):
        for i, _ in MUTATOR:
            self.assertEqual(neomutator_str2int(neomutator_int2str(i)), i)
        self.assertEqual(neomutator_int2str(999), "Unknown Mutator")
        self.assertEqual(neomutator_str2int("invalid"), -1)
        self.assertEqual(neomutator_str2int(None), -1)

    def test_neoselector_mappings(self):
        for i, s in SELECTOR:
            self.assertEqual(neoselector_int2str(i), s)
            self.assertEqual(neoselector_str2int(s), i)

    def test_neoselector_roundtrip_and_unknowns(self):
        for i, _ in SELECTOR:
            self.assertEqual(neoselector_str2int(neoselector_int2str(i)), i)
        self.assertEqual(neoselector_int2str(999), "Unknown Selector")
        self.assertEqual(neoselector_str2int("bad"), -1)
        self.assertEqual(neoselector_str2int(None), -1)

    def test_neosolver_mappings(self):
        for i, s in SOLVER:
            self.assertEqual(neosolver_int2str(i), s)
            self.assertEqual(neosolver_str2int(s), i)

    def test_neosolver_roundtrip_and_unknowns(self):
        for i, _ in SOLVER:
            self.assertEqual(neosolver_str2int(neosolver_int2str(i)), i)
        self.assertEqual(neosolver_int2str(999), "Unknown Solver")
        self.assertEqual(neosolver_str2int("unknown"), -1)
        self.assertEqual(neosolver_str2int(None), -1)


if __name__ == "__main__":
    unittest.main()