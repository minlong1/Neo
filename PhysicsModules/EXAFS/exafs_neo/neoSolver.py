import numpy as np

from PhysicsModules.EXAFS.exafs_neo.exafs_pop import NeoPopulations
from PhysicsModules.EXAFS.exafs_neo.neoPars import NeoPars
from PhysicsModules.EXAFS.exafs_neo.utils import NeoLogger


class NeoSolverBase:
    def __init__(self, exafs_pars: NeoPars, logger: NeoLogger):
        """
        Initialize the selector base class
        :param exafs_pars:
        :param logger:
        """
        self.logger = logger
        self.exafs_pars = exafs_pars

        self.sol_list = []

    def solve(self, pops, selector, crossover, mutator, exafs_pars):
        pass

    def __str__(self):
        # return f"Top Percentage: {100 * self.nBest_Percent}%, Lucky: {100 * self.nLucky_Percent}%"
        return "Neo Solver"


class NeoSolver_GA(NeoSolverBase):
    """
    Standard GA algorithm solver
    """

    def __init__(self, exafs_pars, logger):
        super().__init__(exafs_pars, logger)
        self.solver_type = 0
        self.solver_operator = "Genetic Algorithm"

    def solve(self, pops, selector, crossover, mutator, exafs_pars):
        selector.select(pops)
        crossover.crossover(pops)
        mutator.mutate(pops)
        pops.eval_population()


class NeoSolver_GA_Rechenberg(NeoSolverBase):
    """
    Standard GA with Rechenberg addition
    """

    def __init__(self, exafs_pars, logger):
        super().__init__(exafs_pars, logger)
        self.solver_type = 1
        self.solver_operator = "Genetic Algorithm with Rechenberg"
        self.diffCounter = 0

    def solve(self, pops, selector, crossover, mutator, neo_pars):
        selector.select(pops)
        crossover.crossover(pops)
        self.rechenberg_mutation(neo_pars)
        mutator.mutate(pops)
        pops.eval_population()

    def rechenberg_mutation(self, neo_pars: NeoPars):
        # Recehenberg mutation
        diffCounter = neo_pars.runPars.diffCounter
        bestDiff = np.abs(
            neo_pars.bestFitPars.globBestVal - neo_pars.bestFitPars.currBestVal
        )
        if neo_pars.runPars.currGen > 20:
            if bestDiff < 0.1:
                diffCounter += 1
            else:
                diffCounter -= 1
            # Clip between 0 and None
            diffCounter = np.clip(diffCounter, 0, None)

            if (abs(diffCounter) / float(neo_pars.runPars.currGen)) > 0.2:
                if (neo_pars.mutPars.mutChance + 0.0025) < 1.0:
                    neo_pars.mutPars.mutChance += 0.0025
                    neo_pars.mutPars.mutChance = abs(neo_pars.mutPars.mutChance)
            elif (abs(diffCounter) / float(neo_pars.runPars.currGen)) < 0.2:
                if (neo_pars.mutPars.mutChance - 0.0025) > 0:
                    neo_pars.mutPars.mutChance -= 0.0025
                    # neo_pars.mutPars.mutChance = abs(neo_pars.mutPars.mutChance)
            # Clip between 0 and 100%
            neo_pars.mutPars.mutChance = np.clip(neo_pars.mutPars.mutChance, 0, 1.0)

        neo_pars.runPars.diffCounter = diffCounter
        return neo_pars


class NeoSolver_DE(NeoSolverBase):
    """
    Standard Differential Evolution
    """

    def __init__(self, exafs_pars, logger):
        super().__init__(exafs_pars, logger)
        self.solver_type = 2
        self.solver_operator = "Differential Evolution"

    def solve(self, pops, selector, crossover, mutator, exafs_pars):
        pass


class NeoSolver:
    def __init__(self, logger=None):
        """
        Neo Selector
        :param NeoLogger logger: logger for Neo
        """
        self.solver_operator = None
        self.logger = logger
        self.solver_type = None
        self.exafs_pars = None

    def initialize(self, exafs_pars):
        """
        Initialize the Selector
        :param exafs_pars:
        :return:
        """
        self.exafs_pars = exafs_pars
        # self.solver_type = exafs_pars.selPars.selOpt
        self.solver_type = exafs_pars.solPars.solOpt
        if self.solver_type == 0:
            self.solver_operator = NeoSolver_GA(exafs_pars, logger=self.logger)
        elif self.solver_type == 1:
            self.solver_operator = NeoSolver_GA_Rechenberg(
                exafs_pars, logger=self.logger
            )
        elif self.solver_type == 2:
            self.solver_operator = NeoSolver_DE(exafs_pars, logger=self.logger)
        else:
            self.solver_operator = NeoSolverBase(exafs_pars, logger=self.logger)
            raise ValueError("Invalid selector type, returning standard selector type.")

    def solve(self, pops, selector, crossover, mutator, exafs_pars):
        """
        Perform the actual selection
        :param exafs_pars:
        :param selector:
        :param mutator:
        :param crossover:
        :param NeoPopulation pops:
        :return:
        """
        if self.solver_operator is None:
            raise ValueError("Solver is not initialized")
        else:
            return self.solver_operator.solve(
                pops, selector, crossover, mutator, exafs_pars
            )

    def __str__(self):
        if self.solver_operator is None:
            return "None Mutator selected"
        else:
            return f"Selector Type: {self.solver_type}, {self.solver_operator}"


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
        "solver_type": 1,
    }
    exafs_NeoPars = NeoPars()
    exafs_NeoPars.read_inputs(inputs_pars)

    neo_population = NeoPopulations(exafs_NeoPars)
    neo_population.initialize_populations()

    exafs_solver = NeoSolver()
    exafs_solver.initialize(exafs_pars=exafs_NeoPars)
    # exafs_selector.solve(neo_population)
    print(exafs_solver)
