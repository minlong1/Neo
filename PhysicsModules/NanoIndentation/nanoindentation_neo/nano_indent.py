"""
Nano Neo runner: wires the NanoIndentationProblem to the shared NEO solvers.

Mirrors the EXAFS module pattern (read -> setup -> run); the GA itself lives
in Solvers. Historical .ini `mutated_options` values map onto Solvers
operators as:

    0 -> Rechenberg-adaptive GA, mutate per individual   (original default)
    1 -> GA, mutate per gene ("per path")
    2 -> GA, Metropolis mutation
    3 -> GA, bounded mutation
"""

import csv
import os
import pathlib
import time
from typing import Dict, Union

import numpy as np

from Solvers import get_solver
from Solvers.core import SolverResult

from PhysicsModules.NanoIndentation.nanoindentation_neo.helper import (
    NanoLogger,
    banner,
    bcolors,
)
from PhysicsModules.NanoIndentation.nanoindentation_neo.parser import parse_input_file
from PhysicsModules.NanoIndentation.nanoindentation_neo.problem import (
    NanoIndentationProblem,
)

# (solver key, Solvers mutOpt) per historical mutated_options value
_MUT_OPT_MAP = {
    0: ("GA_RECHENBERG", 0),
    1: ("GA", 1),
    2: ("GA", 3),
    3: ("GA", 4),
}


class NanoNeo:
    """Top-level Nano Neo run: read -> setup -> run, mirroring
    ExafsNeo/XPS_GA's lifecycle. `nano_read` loads the .ini or a parameter
    dict, `nano_setup` builds the NanoIndentationProblem and solver,
    `run` drives the generation loop and writes per-generation outputs."""

    def __init__(self, verbose_lvl: int = 5):
        self.logger = NanoLogger()
        self.verbose_lvl = verbose_lvl
        self.input_parameters = None
        self.problem = None
        self.solver = None
        self.output_path = None
        self.data_path = None

    def nano_read(self, filepath: Union[str, pathlib.Path] = None, input_parameters: Dict = None) -> None:
        """
        Read the input file (or a parameter dict) into Nano Neo parameters
        """
        if filepath is not None:
            file_path = os.path.join(os.getcwd(), str(filepath))
            self.input_parameters = parse_input_file(file_path)

        if input_parameters is not None:
            self.input_parameters = input_parameters

        if self.input_parameters is None:
            raise ValueError("No input parameters are given")

    def nano_setup(self) -> None:
        """
        Setup the problem, output files, and solver
        """
        pars = self.input_parameters

        self._initialize_outputs(pars["output_file"])
        self.logger.print(banner())

        self.problem = NanoIndentationProblem(
            data_file=pars["data_file"],
            data_cutoff=pars.get("data_cutoff", (0.1, 0.9)),
            npaths=pars.get("npaths", 1),
            A_range=pars.get("A_range", ()),
            hf_range=pars.get("hf_range", ()),
            m_range=pars.get("m_range", ()),
            fits=pars.get("fits", "OliverPharr"),
        )

        solver_key, mut_opt = _MUT_OPT_MAP[pars.get("mutOpt", 0)]
        solver_cls = get_solver(solver_key)
        self.solver = solver_cls(
            self.problem,
            options={
                "nPops": pars.get("nPops", 100),
                "nGen": pars.get("nGen", 100),
                "selOpt": 0,
                "croOpt": 0,  # uniform crossover, as in nano_indent
                "mutOpt": mut_opt,
                "mutChance": pars.get("mutChance", 0.2),
                "nBestSample": pars.get("nBestSample", 0.2),
                "nLuckSample": pars.get("nLuckSample", 0.1),
            },
            logger=self.logger,
        )

    def run(self) -> SolverResult:
        """
        Run the fit; owns the generation loop for per-generation
        logging/output (same pattern as the EXAFS module).
        """
        solver = self.solver
        self._run_verbose_start()

        solver.initialize()
        for _ in range(solver.n_gen):
            solver.state.start_gen()
            gen_start = time.time()
            solver.step()
            solver.result.collect(solver.state, solver.population)
            self._gen_time = time.time() - gen_start
            self._run_verbose_gen(solver)
            self._output_generations(solver)
            solver.state.end_gen()
        solver.problem.on_run_end(solver.state, solver.population)

        self._run_verbose_end(solver)
        if self.input_parameters.get("printGraph", False):
            self.plot(solver.state.globBestInd)
        return solver.result

    def plot(self, individual):
        import matplotlib.pyplot as plt

        problem = self.problem
        plt.figure()
        raw = problem.data_obj.get_raw_data()
        plt.plot(raw[:, 0], raw[:, 1], "b-.", label="raw")
        plt.plot(problem.x_slice, problem.y_slice, "o--", label="data")
        plt.plot(problem.x_slice, problem.model_total(individual.genes), "r--", label="model")
        plt.legend()
        plt.show()

    # ------------------------- internals -------------------------

    def _initialize_outputs(self, output_file: str) -> None:
        self.output_path = os.path.join(os.getcwd(), output_file)
        self._check_if_exists(self.output_path)
        with open(self.output_path, "a+") as f:
            f.write("Gen,TPS,FITTNESS,CURRFIT,CURRIND,BESTFIT,BESTIND\n")

        self.data_path = os.path.splitext(self.output_path)[0] + "_data.csv"
        self._check_if_exists(self.data_path)

        log_path = os.path.splitext(self.output_path)[0] + ".log"
        self._check_if_exists(log_path)
        self.logger.initialize_logging(log_path)

    @staticmethod
    def _check_if_exists(path_file) -> None:
        if os.path.exists(path_file):
            os.remove(path_file)
        pathlib.Path(path_file).parent.mkdir(parents=True, exist_ok=True)

    def _run_verbose_start(self) -> None:
        log = self.logger.print
        log("-----------Inputs File Stats---------------")
        log(f"{bcolors.BOLD}File{bcolors.ENDC}: {self.problem.data_file}")
        log(f"{bcolors.BOLD}File Type{bcolors.ENDC}: {self.problem.data_obj._ftype}")
        log(f"{bcolors.BOLD}Output{bcolors.ENDC}: {self.output_path}")
        log(f"{bcolors.BOLD}Population{bcolors.ENDC}: {self.solver.n_pops}")
        log(f"{bcolors.BOLD}Num Gen{bcolors.ENDC}: {self.solver.n_gen}")
        log(f"{bcolors.BOLD}Solver{bcolors.ENDC}: {type(self.solver).name}")
        log("-------------------------------------------")

    def _run_verbose_gen(self, solver) -> None:
        if self.verbose_lvl < 5:
            return
        state = solver.state
        log = self.logger.print
        log("---------------------------------------------------------")
        log(f"{bcolors.BOLD}Gen: {bcolors.ENDC}{state.currGen}")
        with np.printoptions(precision=5, suppress=True):
            log(f"{bcolors.BOLD}Best fit: {bcolors.OKBLUE}{state.currBestVal}{bcolors.ENDC}")
            log("Best fit combination:\n" + str(np.asarray(state.currBestInd.genes)))
            log(f"{bcolors.BOLD}History Best: {bcolors.OKBLUE}{state.globBestVal}{bcolors.ENDC}")
        log(f"Mutation Chance: {solver.mutator.mut_chance:.4f}")
        log("Time: " + str(round(self._gen_time, 5)) + "s")

    def _run_verbose_end(self, solver) -> None:
        log = self.logger.print
        log("-----------Output Stats---------------")
        log(f"{bcolors.BOLD}Total Time(s){bcolors.ENDC}: {round(solver.state.tt, 4)}")
        log(f"{bcolors.BOLD}Final Fitness{bcolors.ENDC}: {solver.state.globBestVal}")
        log("-------------------------------------------")

    def _output_generations(self, solver) -> None:
        """
        Output generation results into the run file and the best-fit data file
        """
        state = solver.state
        with open(self.output_path, "a") as f1:
            f1.write(
                f"{state.currGen},{self._gen_time},"
                f"{state.currBestVal},{list(state.currBestInd.genes)},"
                f"{state.globBestVal},{list(state.globBestInd.genes)}\n"
            )
        with open(self.data_path, "a") as f2:
            writer = csv.writer(f2)
            best = np.asarray(state.globBestInd.genes).reshape(self.problem.npaths, 3)
            for row in best:
                writer.writerow((row[0], row[1], row[2]))
            f2.write("#################################\n")
