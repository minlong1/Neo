import copy
import numpy as np
from PhysicsModules.EXAFS.exafs_neo.exafs_pop import NeoPopulations, fitness
from PhysicsModules.EXAFS.exafs_neo.neoPars import NeoPars


class ExafsMutatorBase:
    def __init__(self, exafs_pars, logger):
        self.logger = logger
        self.exafs_pars = exafs_pars
        self.mutOpt = self.exafs_pars.mutPars.mutOpt
        self.mutChance = self.exafs_pars.mutPars.mutChance
        self.mutChanceE0 = self.exafs_pars.mutPars.mutChanceE0
        self.mutType = None

    def mutate(self, pops):
        pass

    def __str__(self):
        return f"mutation chance: {self.mutChance}%, mutation chance E0: {self.mutChanceE0}%"


class ExafsMutator_PerIndividual(ExafsMutatorBase):
    def __init__(self, exafs_pars, logger):
        super().__init__(exafs_pars, logger)
        self.mutOpt = 1

        self.mutType = "Mutate Per Individual"

    def mutate(self, pops):
        for i, _ in enumerate(pops.population):
            if np.random.random() < self.mutChance:
                new_ind = pops.generate_individual()
                pops.population[i] = new_ind
                self.exafs_pars.mutPars.nmut += 1


class ExafsMutator_PerPath(ExafsMutatorBase):
    def __init__(self, exafs_pars, logger):
        super().__init__(exafs_pars, logger)
        self.mutOpt = 2

        self.mutType = "Mutate Per Path"

    def mutate(self, pops):
        for i, individual in enumerate(pops.population):
            if np.random.random() < self.mutChance:
                pops.population[i].mutate_paths(self.mutChance)
                self.exafs_pars.mutPars.nmut += 1


class ExafsMutator_PerTrait(ExafsMutatorBase):
    def __init__(self, exafs_pars, logger):
        super().__init__(exafs_pars, logger)
        self.mutOpt = 3

        self.mutType = "Mutate Per Trait"

    def mutate(self, pops):
        for i, pop in enumerate(pops.population):
            for j, trait in enumerate(pop):
                if np.random.random() < self.mutChance:
                    pass

                    # pops.population[i].mutate_paths(self.mutChance)
                    # self.exafs_pars.mutPars.nmut += 1


class ExafsMutator_Metropolis(ExafsMutatorBase):
    def __init__(self, exafs_pars, logger):
        super().__init__(exafs_pars, logger)
        self.mutOpt = 4

        self.mutType = "Mutate Metropolis"

    def mutate(self, pops):
        for i, indi in enumerate(pops.population):
            if np.random.random() < self.mutChance:
                nmutate_success = 0

                og_indi = copy.deepcopy(indi)
                og_score = fitness(self.exafs_pars, og_indi)
                mut_indi = copy.deepcopy(indi)
                mut_indi.mutate_paths(self.mutChance)
                mut_score = fitness(self.exafs_pars, mut_indi)
                # T = - self.bestDiff / np.log(1 - (self.genNum / self.ngen))
                T = -self.exafs_pars.bestFitPars.bestDiff / np.log(
                    1
                    - (self.exafs_pars.runPars.currGen / self.exafs_pars.fixedPars.nGen)
                )
                if mut_score < og_score:
                    nmutate_success += 1
                    newIndi = mut_indi
                elif np.exp(-(mut_score - og_score) / T) > np.random.uniform():
                    nmutate_success += 1
                    newIndi = mut_indi
                else:
                    newIndi = og_indi

                pops.population[i] = newIndi


class ExafsMutator_Bounded(ExafsMutatorBase):
    # TODO: This is not working...
    def __init__(self, exafs_pars, logger):
        super().__init__(exafs_pars, logger)
        self.mutOpt = 5

        self.mutType = "Mutate Bounded Per Range"

    def mutate(self, pops):
        def delta_fun(t, delta_val):
            rnd = np.random.random()
            return delta_val * (
                1 - rnd ** (1 - (t / self.exafs_pars.fixedPars.nGen)) ** 5
            )

        for i, indi in enumerate(pops.population):
            if np.random.random() < self.mutChance:
                og_indi = copy.deepcopy(indi)
                og_data = og_indi.get_var()
                for j, path in enumerate(og_data):
                    print(j, path)
                    arr = np.random.randint(2, size=3)
                    for k in range(len(arr)):
                        new_path = []
                        val = path[k]
                        if arr[k] == 0:
                            UP = self.exafs_pars.exafsRangePars.pathrange_pars[
                                j
                            ].get_lim()[k + 1][1]
                            del_val = delta_fun(
                                self.exafs_pars.runPars.currGen, UP - val
                            )
                            val = val + del_val
                        if arr[k] == 1:
                            LB = self.exafs_pars.exafsRangePars.pathrange_pars[
                                j
                            ].get_lim()[k + 1][0]
                            del_val = delta_fun(
                                self.exafs_pars.runPars.currGen, val - LB
                            )
                            val = val - del_val
                        new_path.append(val)
                    indi.set_path(j, new_path[0], new_path[1], new_path[2])


class NeoMutator:
    def __init__(self, logger=None):
        self.mutator = None
        self.logger = logger
        self.mutator_type = None
        self.exafs_pars = None
        self.mutator_score = 0  # TODO: need to check if this is needed

    def initialize(self, exafs_pars):
        self.exafs_pars = exafs_pars

        self.mutator_type = exafs_pars.mutPars.mutOpt
        if self.mutator_type == 0:
            self.mutator = ExafsMutator_PerIndividual(
                self.exafs_pars, logger=self.logger
            )
        elif self.mutator_type == 1:
            self.mutator = ExafsMutator_PerPath(self.exafs_pars, logger=self.logger)
        elif self.mutator_type == 2:
            self.mutator = ExafsMutator_PerTrait(self.exafs_pars, logger=self.logger)
        elif self.mutator_type == 3:
            self.mutator = ExafsMutator_Metropolis(self.exafs_pars, logger=self.logger)
        elif self.mutator_type == 4:
            self.mutator = ExafsMutator_Bounded(self.exafs_pars, logger=self.logger)
        elif self.mutator_type == 5:
            self.mutator = ExafsMutator_DE(self.exafs_pars, logger=self.logger)
        else:
            self.mutator = ExafsMutatorBase(self.exafs_pars, logger=self.logger)
            raise ValueError("Invalid mutator type")

        return self.mutator

    def __str__(self):
        if self.mutator is None:
            return "None Mutator selected"
        else:
            return f"Mutator Type: {self.mutator.mutType}, {self.mutator}"

    def mutate(self, pops):
        if self.mutator is None:
            raise ValueError("Mutator is not initialized")
        else:
            if self.exafs_pars.runPars.secondHalf:
                self.mutate_e0(pops)
            self.mutator.mutate(pops)

    def mutate_e0(self, pops):
        """
        Mutate the e0 value in the second half
        """
        if np.random.random() * 100 < self.mutator.mutChanceE0:
            e0 = np.random.choice(self.exafs_pars.exafsRangePars.rangeE0)
            if self.exafs_pars.verbose_lvl >= 5:
                self.logger.print(f"Mutate e0 to: {e0:.3f}")
            for individual in pops.population:
                individual.set_e0(e0)


class ExafsMutator_DE(ExafsMutatorBase):
    """
    Need to rewrite this in the future for separate DE parameters
    """

    def __init__(self, exafs_pars, logger):
        super().__init__(exafs_pars, logger)
        self.mutOpt = 6

        self.mutType = "Mutate DE"

    def mutate(self, pops) -> list:
        pass


if __name__ == "__main__":
    inputs_pars = {
        "data_file": "../path_files/Cu/cu_10k.xmu",
        "output_file": "",
        "feff_file": "../path_files/Cu/path_75/feff",
        "kmin": 0.95,
        "kmax": 9.775,
        "kweight": 3.0,
        "pathrange": [1, 2, 3, 4, 5],
        "deltak": 0.05,
        "rbkg": 1.1,
        "bkgkw": 1.0,
        "bkgkmax": 15.0,
        "mut_options": 1,
    }
    exafs_NeoPars = NeoPars()
    exafs_NeoPars.read_inputs(inputs_pars)
    # Initialize the population
    neo_population = NeoPopulations(exafs_NeoPars)
    neo_population.initialize_populations()
    print(neo_population.population[0].get_var())

    # Initialize the mutator
    exafs_mutator = NeoMutator()
    exafs_mutator.initialize(exafs_pars=exafs_NeoPars)
    print(exafs_mutator)
    exafs_mutator.mutate(neo_population)
    print(neo_population.population[0].get_var())
