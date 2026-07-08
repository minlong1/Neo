import copy
import operator
from typing import List, Any

import numpy as np
from attrs import define, field
from larch.xafs import feffdat

from PhysicsModules.EXAFS.exafs_neo.individual import Individual
from PhysicsModules.EXAFS.exafs_neo.neoPars import NeoPars
from PhysicsModules.EXAFS.exafs_neo.problem import EXAFSProblem, NeoRunStateView


def fitness(exafs_neo_pars, ind_obj, return_tot=False):
    """
    Evaluate fitness of an individual using vectorized operations.
    """
    y_total = np.zeros(401)

    # Cast intervalK to an integer NumPy array once for fast indexing
    idx = np.array(exafs_neo_pars.exafsPars.intervalK, dtype=int)

    larch = exafs_neo_pars.exafsPathPars.mylarch
    kweight = exafs_neo_pars.exafsPars.kweight

    # Fetch e0 once if it is constant across all paths for this individual
    e0_val = ind_obj.get_e0()

    for i in range(exafs_neo_pars.exafsPars.npath):
        pathname = exafs_neo_pars.exafsPathPars.pathname[i]
        path = exafs_neo_pars.exafsPathPars.pathDictionary.get(pathname)

        # Call get_path(i) once to save repeated method calls
        path_params = ind_obj.get_path(i)

        path.e0 = e0_val
        path.s02 = path_params[0]
        path.sigma2 = path_params[2]
        path.deltar = path_params[3]

        feffdat.path2chi(path, _larch=larch)

        # Vectorized addition: Update only the target indices at once
        y_total[idx] += path.chi[idx]

    # Vectorized loss computation
    g_k_idx = exafs_neo_pars.exafsPathPars.g.k[idx]
    exp_idx = exafs_neo_pars.exafsPathPars.exp[idx]

    # Calculate the weighted difference, square it, and sum it across all indices
    weighted_diff = (y_total[idx] - exp_idx) * (g_k_idx ** kweight)
    loss = np.sum(weighted_diff ** 2)

    if return_tot:
        return loss, y_total

    return loss


@define
class NeoPopulations:
    """EXAFS's population container, scored via the module-level `fitness`
    (larch path2chi + k-weighted residual). Exposes `.problem`/`.state` as
    Solvers-facing views so generic `Solvers.ga`/`Solvers.de` operators can
    run against it without knowing about larch or EXAFS paths."""

    exafs_NeoPars: NeoPars = None
    population: List = field(factory=list)
    population_sorted: List = field(factory=list)
    population_score: List = field(factory=list)
    population_perf: dict = field(factory=dict)
    next_population: List = field(factory=list)
    # Solvers-facing views, so generic operators can run against this
    # population (pops.problem.fitness, pops.state.currGen, ...)
    problem: EXAFSProblem = None
    state: NeoRunStateView = None

    def __attrs_post_init__(self):
        if self.exafs_NeoPars is not None:
            self.problem = EXAFSProblem(self.exafs_NeoPars)
            self.state = NeoRunStateView(self.exafs_NeoPars)

    def generate_individual(self) -> Individual:
        if not self.exafs_NeoPars.runPars.secondHalf:
            e0 = np.random.choice(self.exafs_NeoPars.exafsRangePars.rangeE0)
        else:
            e0 = self.exafs_NeoPars.bestFitPars.bestE0
        return Individual(
            self.exafs_NeoPars.exafsPars.npath,
            self.exafs_NeoPars.exafsPathPars.pathDictionary,
            self.exafs_NeoPars.exafsRangePars.pathrange_pars,
            self.exafs_NeoPars.exafsPathPars.path_lists,
            e0,
            self.exafs_NeoPars.exafsPathPars.pathname,
        )

    def eval_population(self, replace=True, sorting=True) -> List:
        score = []
        population_perf = {}

        for i, individual in enumerate(self.population):
            temp_score = fitness(self.exafs_NeoPars, individual)
            score.append(temp_score)

            population_perf[individual] = temp_score
        if sorting:
            self.population_sorted = sorted(
                population_perf.items(), key=operator.itemgetter(1), reverse=False
            )
        if replace:
            # self.currBestFit = list(self.population_sorted[0])
            self.__replace_bestfit()
        return score

    def initialize_populations(self) -> None:
        """
        Initialize populations
        :return:
        """
        for i in range(self.exafs_NeoPars.fixedPars.nPops):
            self.population.append(self.generate_individual())

        self.eval_population()

    def __getitem__(self, item) -> list[Any] | Any:
        return self.population_sorted[item]

    def __replace_bestfit(self) -> None:
        self.exafs_NeoPars.bestFitPars.currBestInd = self.population_sorted[0][0]
        self.exafs_NeoPars.bestFitPars.currBestVal = self.population_sorted[0][1]

        # if this is positive, it means we have a better solution
        delta = (
            self.exafs_NeoPars.bestFitPars.globBestVal
            - self.exafs_NeoPars.bestFitPars.currBestVal
        )
        # Check of cutoff minima improvement.
        if np.abs(delta) > 0.01:
            self.exafs_NeoPars.bestFitPars.bestDiff = delta
        else:
            self.exafs_NeoPars.bestFitPars.bestDiff = 0

        if (
            self.exafs_NeoPars.bestFitPars.currBestVal
            < self.exafs_NeoPars.bestFitPars.globBestVal
        ):
            self.exafs_NeoPars.bestFitPars.globBestInd = (
                self.exafs_NeoPars.bestFitPars.currBestInd
            )
            self.exafs_NeoPars.bestFitPars.globBestVal = (
                self.exafs_NeoPars.bestFitPars.currBestVal
            )

    def optimize_e0(self) -> None:
        # TODO: Revisit this
        #  if mess == None:
        #      self.logger.info(
        #          "Finished First Half of Generation, Optimizing E0...")
        #  else:
        #      self.logger.info(mess)

        curr_ind = copy.deepcopy(self.exafs_NeoPars.bestFitPars.globBestInd)
        curr_score = copy.deepcopy(self.exafs_NeoPars.bestFitPars.globBestVal)
        curr_e0 = curr_ind.get_e0()
        for i in self.exafs_NeoPars.exafsRangePars.rangeE0_large:
            curr_ind.set_e0(i)
            fit_score = fitness(self.exafs_NeoPars, curr_ind)
            if fit_score < curr_score:
                curr_e0 = i
                curr_score = fit_score
            # listOfX.append(i)
            # listOfY.append(fit)
        # self.logger.info("Continue With E0= " + str(round(curr_e0, 3)))
        new_e0 = curr_e0
        self.exafs_NeoPars.bestFitPars.bestE0 = new_e0
        # TODO: revisit this!
        #  Reset Mutation Chance??
        #  self.mut_chance_e0 = 0
        self.exafs_NeoPars.bestFitPars.globBestInd.set_e0(new_e0)
        self.exafs_NeoPars.bestFitPars.globBestVal = curr_score

        for i in self.population:
            i.set_e0(new_e0)

        self.exafs_NeoPars.bestFitPars.bestE0 = new_e0


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
    neo_population.eval_population()
