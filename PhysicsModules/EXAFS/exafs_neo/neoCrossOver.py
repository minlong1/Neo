import numpy as np

from PhysicsModules.EXAFS.exafs_neo.exafs_pop import NeoPopulations
from PhysicsModules.EXAFS.exafs_neo.neoPars import NeoPars


class EXAFS_CrossoverBase:
    def __init__(self, exafs_pars, logger=None):
        self.logger = logger
        self.exafs_pars = exafs_pars
        self.croOpt = self.exafs_pars.crossPars.croOpt
        self.croType = None

    def crossover(self, pops, individual1, individual2):
        pass

    def __str__(self):
        return f"Crossover Option: {self.croType}"


class EXAFS_UniformCrossover(EXAFS_CrossoverBase):
    def __init__(self, exafs_pars, logger):
        super().__init__(exafs_pars, logger)
        self.croOpt = 0
        self.croType = "Uniform Crossover"

    def crossover(self, pops, individual1, individual2):
        child = pops.generate_individual()
        if np.random.randint(0, 1):
            child.set_e0(individual1.get_e0())
        else:
            child.set_e0(individual2.get_e0())

        for i in range(self.exafs_pars.exafsPathPars.npath):
            individual1_path = individual1.get_path(i)
            individual2_path = individual2.get_path(i)

            temp_path = []
            for path_pars in range(4):
                if np.random.randint(0, 2):
                    temp_path.append(individual1_path[path_pars])
                else:
                    temp_path.append(individual2_path[path_pars])

            child.set_path(i, temp_path[0], temp_path[2], temp_path[3])

        return child


class EXAFS_SinglePointCrossover(EXAFS_CrossoverBase):
    def __init__(self, exafs_pars, logger):
        super().__init__(exafs_pars, logger)
        self.croOpt = 1
        self.croType = "Single Point Crossover"

    def crossover(self, pops, individual1, individual2, co_point=1):
        # prevent overflow by clipping to 0 and 3 for the path
        co_point = np.clip(co_point, 0, 3)
        child = pops.generate_individual()

        if np.random.randint(0, 1):
            child.set_e0(individual1.get_e0())
        else:
            child.set_e0(individual2.get_e0())

        for i in range(self.exafs_pars.exafsPathPars.npath):
            individual1_path = individual1.get_path(i)
            individual2_path = individual2.get_path(i)

            temp_path = []
            for j in range(4):
                if j < co_point:
                    temp_path.append(individual1_path[j])
                else:
                    temp_path.append(individual2_path[j])

            child.set_path(i, temp_path[0], temp_path[2], temp_path[3])

        return child


class EXAFS_DualPointCrossover(EXAFS_CrossoverBase):
    def __init__(self, exafs_pars, logger):
        super().__init__(exafs_pars, logger)
        self.croOpt = 2
        self.croType = "Dual Point Crossover"

    def crossover(self, pops, individual1, individual2):
        pass


class EXAFS_ArithmeticCrossover(EXAFS_CrossoverBase):
    def __init__(self, exafs_pars, logger):
        super().__init__(exafs_pars, logger)
        self.croOpt = 3
        self.croType = "Arithmetic Crossover"

    def crossover(self, pops, individual1, individual2):
        child = pops.generate_individual()
        if np.random.randint(0, 1):
            child.set_e0(individual1.get_e0())
        else:
            child.set_e0(individual2.get_e0())

        for i in range(self.exafs_pars.exafsPathPars.npath):
            individual1_path = individual1.get_path(i)
            individual2_path = individual2.get_path(i)

            temp_path = []
            for j in range(4):
                ind_1 = np.random.randint(0, 2)
                ind_2 = np.random.randint(0, 2)
                if np.logical_and(ind_1, ind_2):
                    temp_path.append(individual1_path[j])
                else:
                    temp_path.append(individual2_path[j])

            child.set_path(i, temp_path[0], temp_path[2], temp_path[3])

        return child


class EXAFS_OrCrossover(EXAFS_CrossoverBase):
    def __init__(self, exafs_pars, logger):
        super().__init__(exafs_pars, logger)
        self.croOpt = 4
        self.croType = "Or Crossover"

    def crossover(self, pops, individual1, individual2):
        child = pops.generate_individual()
        if np.random.randint(0, 1):
            child.set_e0(individual1.get_e0())
        else:
            child.set_e0(individual2.get_e0())

        for i in range(self.exafs_pars.exafsPathPars.npath):
            individual1_path = individual1.get_path(i)
            individual2_path = individual2.get_path(i)

            temp_path = []
            for j in range(4):
                ind_1 = np.random.randint(0, 2)
                ind_2 = np.random.randint(0, 2)
                if np.logical_or(ind_1, ind_2):
                    temp_path.append(individual1_path[j])
                else:
                    temp_path.append(individual2_path[j])

            child.set_path(i, temp_path[0], temp_path[2], temp_path[3])

        return child


class EXAFS_AverageCrossOver(EXAFS_CrossoverBase):
    def __init__(self, exafs_pars, logger):
        super().__init__(exafs_pars, logger)
        self.croOpt = 5
        self.croType = "Average Crossover"

    def crossover(self, pops, individual1, individual2):
        child = pops.generate_individual()
        if np.random.randint(0, 1):
            child.set_e0(individual1.get_e0())
        else:
            child.set_e0(individual2.get_e0())

        for i in range(self.exafs_pars.exafsPathPars.npath):
            individual1_path = individual1.get_path(i)
            individual2_path = individual2.get_path(i)

            temp_path = []
            for j in range(4):
                temp_path.append((individual1_path[j] + individual2_path[j]) / 2)

            # TODO check if this values goes out of bound
            child.set_path(i, temp_path[0], temp_path[2], temp_path[3])

        return child


class NeoCrossover:
    def __init__(self, logger=None):
        self.logger = logger
        self.exafs_pars = None
        self.crossover_type = None
        self.crossover_operator = None
        self.crossover_score = 0  # TODO: maybe implement this?

    def initialize(self, exafs_pars):
        self.exafs_pars = exafs_pars

        self.crossover_type = exafs_pars.crossPars.croOpt
        if self.crossover_type == 0:
            self.crossover_operator = EXAFS_UniformCrossover(
                exafs_pars, logger=self.logger
            )
        elif self.crossover_type == 1:
            self.crossover_operator = EXAFS_SinglePointCrossover(
                exafs_pars, logger=self.logger
            )
        elif self.crossover_type == 2:
            self.crossover_operator = EXAFS_DualPointCrossover(
                exafs_pars, logger=self.logger
            )
        elif self.crossover_type == 3:
            self.crossover_operator = EXAFS_ArithmeticCrossover(
                exafs_pars, logger=self.logger
            )
        elif self.crossover_type == 4:
            self.crossover_operator = EXAFS_OrCrossover(exafs_pars, logger=self.logger)
        elif self.crossover_type == 5:
            self.crossover_operator = EXAFS_AverageCrossOver(
                exafs_pars, logger=self.logger
            )
        else:
            self.crossover_operator = EXAFS_CrossoverBase(
                exafs_pars, logger=self.logger
            )
            raise ValueError(
                "Invalid crossover type, returning standard crossover type."
            )

    def __str__(self):
        if self.crossover_operator is None:
            return "Crossover is not selected"
        else:
            return f"Crossover Type: {self.crossover_type}, {self.crossover_operator}"

    def crossover(self, pops):
        if self.crossover_operator is None:
            raise ValueError("Crossover is not initialized")
        else:
            temp_population = []
            if len(pops.next_population) > 2:
                for _ in range(self.exafs_pars.selPars.nCross):
                    par_ind = np.random.choice(
                        len(pops.next_population), size=2, replace=False
                    )
                    ind1 = pops.next_population[par_ind[0]]
                    ind2 = pops.next_population[par_ind[1]]
                    child = self.crossover_operator.crossover(pops, ind1, ind2)
                    temp_population.append(child)

                pops.next_population.extend(temp_population)
                pops.population = pops.next_population

    def crossover_single(self, pops, ind1, ind2):
        if self.crossover_operator is None:
            raise ValueError("Crossover is not initialized")
        else:
            return self.crossover_operator.crossover(pops, ind1, ind2)


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
        "croOpt": 1,
    }
    exafs_Pars = NeoPars()
    exafs_Pars.read_inputs(inputs_pars)

    neo_population = NeoPopulations(exafs_Pars)
    neo_population.initialize_populations()
    print(neo_population.population[0].get_var())

    crossover_operator = NeoCrossover()
    crossover_operator.initialize(exafs_pars=exafs_Pars)
    crossover_operator.crossover(neo_population)
    print(crossover_operator)
