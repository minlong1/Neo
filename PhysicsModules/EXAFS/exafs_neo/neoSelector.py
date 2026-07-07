"""
EXAFS-facing selector wrapper. The selection algorithms live in
Solvers.ga.selectors; this module adapts them to the historical EXAFS
API (initialize(exafs_pars) / select(pops)).
"""

from Solvers.ga.selectors import create_selector

from PhysicsModules.EXAFS.exafs_neo.utils import NeoLogger


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
        self.selector_operator = create_selector(
            self.selector_type,
            n_pops=exafs_pars.fixedPars.nPops,
            n_best_sample=exafs_pars.selPars.nBestSample,
            n_lucky_sample=exafs_pars.selPars.nLuckSample,
            logger=self.logger,
        )

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
