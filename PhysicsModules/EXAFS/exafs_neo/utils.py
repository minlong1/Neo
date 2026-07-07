import datetime
import logging
import time
from pathlib import Path
from typing import Union, Dict

import numpy as np

from PhysicsModules.common.cli import BaseLogger
from PhysicsModules.common.colors import TermColors
from PhysicsModules.EXAFS.exafs_neo.utils_mapping import (
    neocrossover_int2str,
    neomutator_int2str,
    neoselector_int2str,
    neosolver_int2str,
)


def time_call():
    return time.time()


def raise_error(msg="Error"):
    """
    Raise an error with a message
    :param msg:
    :return:
    """
    raise Exception(msg)


def checkKey(key, dictionary, alt_value=None, logger=None, verbose=False):
    # TODO: Raise checker for alternative response
    if key not in dictionary:
        if verbose:
            warn_str = f"{key} not found! Using default value of {key}: {alt_value}"
            if logger is None:
                print(warn_str)
            else:
                logger.warn(warn_str)
        return alt_value
    else:
        return dictionary[key]


class NeoLogger(BaseLogger):
    def __init__(self):
        super().__init__(name="", level=logging.DEBUG)


def check_if_exists(path_file: Union[str, Path]) -> None:
    """
    Check if the directory exists
    """
    path_file = Path(path_file)
    if path_file.is_file():
        path_file.unlink()
    # Make Directory when its missing
    path_file.parent.mkdir(parents=True, exist_ok=True)


class STRColors(TermColors):
    @staticmethod
    def logger_print_based_on_verbose_lvl(
        logger: NeoLogger, string: str, verbose_lvl: int, des_verbose_lvl: int
    ) -> None:
        if verbose_lvl >= des_verbose_lvl:
            logger.print(string)

    # TODO: Implement verbose_lvl argument that detects if jupyter is running...
    @staticmethod
    def run_verbose_start(logger, exafs_NeoPars, verbose_lvl: int = 5) -> None:
        """
        Visualize the verbose start place
        """
        # logger.print(banner())
        if exafs_NeoPars.fixedPars.debug_mode:
            STRColors.logger_print_based_on_verbose_lvl(
                logger, f"{STRColors.BOLD}DEBUG-MODE{STRColors.ENDC}", verbose_lvl, 5
            )
        STRColors.logger_print_based_on_verbose_lvl(
            logger, "------------Inputs File Stats--------------", verbose_lvl, 5
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Data File{STRColors.ENDC}: {exafs_NeoPars.neoFilePars.data_path}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Feff File{STRColors.ENDC}: {exafs_NeoPars.neoFilePars.front}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Output File{STRColors.ENDC}: {exafs_NeoPars.neoFilePars.output_path}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Log File{STRColors.ENDC}: {exafs_NeoPars.neoFilePars.log_path}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Paths Range File{STRColors.ENDC}: {exafs_NeoPars.neoFilePars.feff_file}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}CSV series{STRColors.ENDC}: {exafs_NeoPars.neoFilePars.multi_data_toggle}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger, "--------------Populations------------------", verbose_lvl, 5
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Population{STRColors.ENDC}: {exafs_NeoPars.fixedPars.nPops}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Num Gen{STRColors.ENDC}: {exafs_NeoPars.fixedPars.nGen}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Best Individuals{STRColors.ENDC}: {exafs_NeoPars.selPars.nBestSample}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Lucky Survivor{STRColors.ENDC}: {exafs_NeoPars.selPars.nLuckSample}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger, "-----------------Paths---------------------", verbose_lvl, 5
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Individual Path{STRColors.ENDC}: {exafs_NeoPars.exafsPathPars.individual_paths}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Num Path{STRColors.ENDC}: {exafs_NeoPars.exafsPathPars.npath}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Paths {STRColors.ENDC}: {exafs_NeoPars.exafsPathPars.path_lists}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Path Optimize{STRColors.ENDC}: {exafs_NeoPars.fixedPars.pathOptimize}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Path Optimize Percent{STRColors.ENDC}: {exafs_NeoPars.fixedPars.pathOptimizePercent}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Path Optimize Only{STRColors.ENDC}: {exafs_NeoPars.fixedPars.pathOptimizeOnly}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger, "-----------------Solvers-------------------", verbose_lvl, 5
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Solver Options{STRColors.ENDC}: {neosolver_int2str(exafs_NeoPars.solPars.solOpt)}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger, "----------------Mutations------------------", verbose_lvl, 5
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Mutations{STRColors.ENDC}: {exafs_NeoPars.mutPars.mutChance}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}E0 Mutations{STRColors.ENDC}: {exafs_NeoPars.mutPars.mutChanceE0}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Mutation Options{STRColors.ENDC}: {neomutator_int2str(exafs_NeoPars.mutPars.mutOpt)}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Selection Options{STRColors.ENDC}: {neoselector_int2str(exafs_NeoPars.selPars.selOpt)}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Crossover Options{STRColors.ENDC}: {neocrossover_int2str(exafs_NeoPars.crossPars.croOpt)}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger, "---------------Larch Paths-----------------", verbose_lvl, 5
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Kmin{STRColors.ENDC}: {exafs_NeoPars.exafsPars.kmin}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Kmax{STRColors.ENDC}: {exafs_NeoPars.exafsPars.kmax}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Kweight{STRColors.ENDC}: {exafs_NeoPars.exafsPars.kweight}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Delta k{STRColors.ENDC}: {exafs_NeoPars.exafsPars.dk}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}RBKG{STRColors.ENDC}: {exafs_NeoPars.exafsPars.rbkg}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}BKG Kw{STRColors.ENDC}: {exafs_NeoPars.exafsPars.bkgkw}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}BKG Kmax{STRColors.ENDC}: {exafs_NeoPars.exafsPars.bkgkmax}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger, "-------------------------------------------", verbose_lvl, 5
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Steady State{STRColors.ENDC}: {exafs_NeoPars.fixedPars.steadyState}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Print Graph{STRColors.ENDC}: {exafs_NeoPars.fixedPars.printGraph}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger, "-------------------------------------------", verbose_lvl, 5
        )

    @staticmethod
    def run_verbose_gen(logger, exafs_NeoPars, neo_population, verbose_lvl=5):
        """
        Verbose generation
        """
        st = time_call()
        currGen = exafs_NeoPars.runPars.currGen
        nGen = exafs_NeoPars.runPars.nGen

        if currGen == nGen // 2 or currGen == nGen:
            STRColors.logger_print_based_on_verbose_lvl(
                logger,
                "---------------------------------------------------------",
                verbose_lvl,
                1,
            )
            STRColors.logger_print_based_on_verbose_lvl(
                logger,
                "Continue With E0= " + str(round(exafs_NeoPars.bestFitPars.bestE0, 3)),
                verbose_lvl,
                5,
            )

        population_sorted = neo_population.population_sorted
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            "---------------------------------------------------------",
            verbose_lvl,
            1,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            datetime.datetime.fromtimestamp(st).strftime("%H:%M:%S")
            + f"{STRColors.BOLD} Gen: {STRColors.ENDC}{exafs_NeoPars.runPars.currGen}",
            verbose_lvl,
            1,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"Best Fit: {STRColors.BOLD}{population_sorted[0][1].round(3)}{STRColors.ENDC}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"2nd Fit: {STRColors.BOLD}{population_sorted[1][1].round(3)}{STRColors.ENDC}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"3rd Fit: {STRColors.BOLD}{population_sorted[2][1].round(3)}{STRColors.ENDC}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"4th Fit: {STRColors.BOLD}{population_sorted[0][1].round(3)}{STRColors.ENDC}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"Last Fit: {STRColors.BOLD}{population_sorted[-1][1].round(3)}{STRColors.ENDC}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"Different from last best fit: {exafs_NeoPars.bestFitPars.bestDiff:.4f}",
            verbose_lvl,
            5,
        )

        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"Glob Best Ind:\n{np.array_str(np.array(exafs_NeoPars.bestFitPars.globBestInd.get()), precision=3)}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"Curr Best Ind:\n{np.array_str(np.array(exafs_NeoPars.bestFitPars.currBestInd.get()), precision=3)}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            STRColors.BOLD
            + "Best fit: "
            + STRColors.OKBLUE
            + str(np.round(exafs_NeoPars.bestFitPars.currBestVal, 3))
            + STRColors.ENDC,
            verbose_lvl,
            5,
        )
        CurrchiR = np.round(
            exafs_NeoPars.bestFitPars.currBestVal
            / (
                len(exafs_NeoPars.exafsPars.intervalK)
                - 3 * exafs_NeoPars.exafsPars.npath
                + 1
            ),
            3,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            STRColors.BOLD
            + "Best fit ChiR: "
            + STRColors.OKBLUE
            + str(CurrchiR)
            + STRColors.ENDC,
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            STRColors.BOLD
            + "History Best: "
            + STRColors.OKBLUE
            + str(np.round(exafs_NeoPars.bestFitPars.globBestVal, 4))
            + STRColors.ENDC,
            verbose_lvl,
            1,
        )
        GlobchiR = exafs_NeoPars.bestFitPars.globBestVal / (
            len(exafs_NeoPars.exafsPars.intervalK)
            - 3 * exafs_NeoPars.exafsPars.npath
            + 1
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            STRColors.BOLD
            + "History Best ChiR: "
            + STRColors.OKBLUE
            + str(np.round(GlobchiR, 4))
            + STRColors.ENDC,
            verbose_lvl,
            1,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            "DiffCounter: " + str(exafs_NeoPars.runPars.diffCounter),
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"Diff %: {(exafs_NeoPars.runPars.diffCounter / exafs_NeoPars.runPars.currGen):.3f}",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"Mutation Chance: {100 * exafs_NeoPars.mutPars.mutChance:.2f}%",
            verbose_lvl,
            5,
        )
        # if exafs_NeoPars.mutPars.mutOpt == 4:
        #     STRColors.logger_print_based_on_verbose_lvl(logger, "Mutation Percentage" + str(np.round(exafs_NeoPars.self.nmutate_success / self.nmutate, 4)),verbose_lvl,5)

        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            "Time: " + str(round(exafs_NeoPars.runPars.currGen_tt, 3)) + "s",
            verbose_lvl,
            5,
        )

    @staticmethod
    def run_verbose_end(logger, exafs_NeoPars, verbose_lvl=5):
        """
        Verbose end
        """

        STRColors.logger_print_based_on_verbose_lvl(
            logger, "---------------Output Stats---------------", verbose_lvl, 1
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Total Time(s){STRColors.ENDC}: {round(exafs_NeoPars.runPars.tt, 4)}",
            verbose_lvl,
            1,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}File{STRColors.ENDC}: {exafs_NeoPars.neoFilePars.data_path}",
            verbose_lvl,
            1,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Path{STRColors.ENDC}: {exafs_NeoPars.exafsPathPars.path_lists}",
            verbose_lvl,
            1,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"{STRColors.BOLD}Final Fittness Score{STRColors.ENDC}: {exafs_NeoPars.bestFitPars.globBestVal:.2f}",
            verbose_lvl,
            1,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger, "-------------------------------------------", verbose_lvl, 1
        )
        if exafs_NeoPars.fixedPars.printGraph:
            # TODO: Reimplement this
            pass
            # self.verbose_graph()

    @staticmethod
    def run_verbose_path_optimize(
        logger, new_optimize_path: Dict, verbose_lvl=5
    ) -> None:
        STRColors.logger_print_based_on_verbose_lvl(
            logger, "---------------Path Optimizations---------------", verbose_lvl, 5
        )
        if not new_optimize_path:
            STRColors.logger_print_based_on_verbose_lvl(
                logger,
                "No path optimizes, consider lowering the cut off percentage.",
                verbose_lvl,
                5,
            )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"""Filter Percentage: {str(round(new_optimize_path["path_optimize_percent"]*100, 2)).zfill(4)} %""",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"""Original Paths: {new_optimize_path["original_paths"]}""",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger,
            f"""Filter Paths: {new_optimize_path["optimize_paths"]}""",
            verbose_lvl,
            5,
        )
        STRColors.logger_print_based_on_verbose_lvl(
            logger, "Optimized Path:", verbose_lvl, 5
        )

        for key, vals in new_optimize_path["optimize_path_lists"].items():
            STRColors.logger_print_based_on_verbose_lvl(
                logger,
                f"\tPath {key}, Contributed {str(round(vals*100, 2)).zfill(4)} %",
                verbose_lvl,
                5,
            )

        STRColors.logger_print_based_on_verbose_lvl(
            logger, "-------------------------------------------", verbose_lvl, 5
        )


class Utils:

    @staticmethod
    def get_sum_path(path_lists, ncomp):
        """Get the sum of each component

        Args:
            path_lists (list): Path lists
            ncomp (int): Number of components

        Returns:
            list: List of sum of each component
        """
        sum_path = []
        if ncomp > 1:
            for i in range(ncomp):
                sum_path.append(len(path_lists[i]))

        return sum_path


if __name__ == "__main__":
    # Testing Logger
    # logger = NeoLogger()
    # logger.initialize_logging('Test.csv')
    # logger.set_loglevel(logging.DEBUG)

    # Test Str Print
    pass
