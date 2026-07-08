"""
EXAFS-facing solver wrapper. The solve step sequencing and the Rechenberg
adaptive-mutation rule live in Solvers; this module adapts them to the
historical EXAFS API (initialize(exafs_pars) / solve(pops, selector, ...)).
"""

from Solvers.de import differential_evolution_step
from Solvers.ga.ga_solver import rechenberg_update

from PhysicsModules.EXAFS.exafs_neo.neoPars import NeoPars
from PhysicsModules.EXAFS.exafs_neo.utils import NeoLogger


class NeoSolverBase:
    def __init__(self, exafs_pars: NeoPars, logger: NeoLogger):
        """
        Initialize the solver base class
        :param exafs_pars:
        :param logger:
        """
        self.logger = logger
        self.exafs_pars = exafs_pars

        self.sol_list = []

    def solve(self, pops, selector, crossover, mutator, exafs_pars):
        pass

    def __str__(self):
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
        # Note: this adjusts mutPars.mutChance, which the active mutator
        # sampled at initialize() time — matching the historical behavior.
        (
            neo_pars.runPars.diffCounter,
            neo_pars.mutPars.mutChance,
        ) = rechenberg_update(
            neo_pars.runPars.currGen,
            neo_pars.runPars.diffCounter,
            neo_pars.bestFitPars.globBestVal,
            neo_pars.bestFitPars.currBestVal,
            neo_pars.mutPars.mutChance,
        )
        return neo_pars


class NeoSolver_DE(NeoSolverBase):
    """
    Differential Evolution. The algorithm itself is
    Solvers.de.differential_evolution_step; NeoPopulations already exposes
    the .problem/.generate_individual()/.eval_population() surface that
    function needs, so no EXAFS-specific DE math lives here. selector,
    crossover, and mutator are accepted for interface parity with the GA
    solvers but unused - DE has its own mutation/crossover built in.
    """

    def __init__(self, exafs_pars, logger):
        super().__init__(exafs_pars, logger)
        self.solver_type = 2
        self.solver_operator = "Differential Evolution"
        self.F = 0.5
        self.CR = 0.9

    def solve(self, pops, selector, crossover, mutator, exafs_pars):
        differential_evolution_step(pops, F=self.F, CR=self.CR)


class NeoSolver:
    """Back-compat shim preserving the historical `NeoSolver().initialize(exafs_pars)`
    / `.solve(pops, selector, crossover, mutator, exafs_pars)` API; picks
    the GA/GA_Rechenberg/DE solver class by `solOpt`."""

    def __init__(self, logger=None):
        """
        Neo Solver
        :param NeoLogger logger: logger for Neo
        """
        self.solver_operator = None
        self.logger = logger
        self.solver_type = None
        self.exafs_pars = None

    def initialize(self, exafs_pars):
        """
        Initialize the Solver
        :param exafs_pars:
        :return:
        """
        self.exafs_pars = exafs_pars
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
            raise ValueError("Invalid solver type, returning standard solver type.")

    def solve(self, pops, selector, crossover, mutator, exafs_pars):
        """
        Perform one generation of the configured solver
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
