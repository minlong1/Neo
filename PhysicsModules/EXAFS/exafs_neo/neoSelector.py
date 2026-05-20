from PhysicsModules.EXAFS.exafs_neo.exafs_pop import NeoPopulations
from PhysicsModules.EXAFS.exafs_neo.neoPars import NeoPars
from PhysicsModules.EXAFS.exafs_neo.utils import NeoLogger


class NeoSelectorBase:
    def __init__(self, exafs_pars: NeoPars, logger: NeoLogger):
        """
        Initialize the selector base class
        :param exafs_pars:
        :param logger:
        """
        self.logger = logger
        self.exafs_pars = exafs_pars
        self.npops = exafs_pars.fixedPars.nPops
        self.nBest_Percent = exafs_pars.selPars.nBestSample
        self.nLucky_Percent = exafs_pars.selPars.nLuckSample

        self.nBest = int(self.nBest_Percent * self.npops)
        self.nLucky = int(self.nLucky_Percent * self.npops)
        self.nStatic = int(self.npops - self.nBest - self.nLucky)

        self.sel_list = []

    def select(self, pops):
        pass

    def __str__(self):
        return f"Top Percentage: {100 * self.nBest_Percent}%, Lucky: {100 * self.nLucky_Percent}%"


class NeoSelector_RouletteWheel(NeoSelectorBase):
    def __init__(self, exafs_pars, logger):
        super().__init__(exafs_pars, logger)
        self.selector_type = 0
        self.selector_operator = "Roulette Wheel"

    def select(self, pops):
        # Create parents
        next_population = []
        for i in range(self.nBest):
            next_population.append(pops.population_sorted[i][0])

        # Create lucky
        for i in range(self.nLucky):
            next_population.append(pops.generate_individual())

        pops.next_population = next_population


class NeoSelector_Tournament(NeoSelectorBase):
    def __init__(self, exafs_pars, logger):
        super().__init__(exafs_pars, logger)
        self.selector_type = 1
        self.selector_operator = "Tournament"

    def select(self, pops):
        pass


class NeoSelector:
    def __init__(self, logger=None):
        """
        Neo Selector
        :param NeoLogger logger: logger for Neo
        """
        self.selector_operator = None
        self.logger = logger
        self.selector_type = None
        self.exafs_pars = None

    def initialize(self, exafs_pars):
        """
        Initialize the Selector
        :param exafs_pars:
        :return:
        """
        self.exafs_pars = exafs_pars
        self.selector_type = exafs_pars.selPars.selOpt
        if self.selector_type == 0:
            self.selector_operator = NeoSelector_RouletteWheel(
                exafs_pars, logger=self.logger
            )
        elif self.selector_type == 1:
            self.selector_operator = NeoSelector_Tournament(
                exafs_pars, logger=self.logger
            )
        else:
            self.selector_operator = NeoSelectorBase(exafs_pars, logger=self.logger)
            raise ValueError("Invalid selector type, returning standard selector type.")

    def select(self, pops):
        """
        Perform the actual selection
        :param NeoPopulation pops:
        :return:
        """
        if self.selector_operator is None:
            raise ValueError("Selector is not initialized")
        else:
            return self.selector_operator.select(pops)

    def __str__(self):
        if self.selector_operator is None:
            return "None Mutator selected"
        else:
            return f"Selector Type: {self.selector_type}, {self.selector_operator}"


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
    }
    exafs_NeoPars = NeoPars()
    exafs_NeoPars.read_inputs(inputs_pars)

    neo_population = NeoPopulations(exafs_NeoPars)
    neo_population.initialize_populations()

    exafs_selector = NeoSelector()
    exafs_selector.initialize(exafs_pars=exafs_NeoPars)
    exafs_selector.select(neo_population)
    print(exafs_selector)
