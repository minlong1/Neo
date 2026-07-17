"""
Astro Neo runner: wires AstroNeoProblem to the shared NEO solvers.

Mirrors the EXAFS/NanoIndentation module pattern (read -> setup -> run); the
GA/DE itself lives in Solvers. The original Astro_Neo demo (astro_neo_test/
AstroNeo.py's __main__) defaulted to solver_type=2 (Differential Evolution),
which is also this module's default -- DE is the better fit for a
15-parameter continuous model like this one.
"""

import csv
import os
import pathlib
import time
from typing import Dict, Union

import numpy as np

from Solvers import get_solver
from Solvers.core import SolverResult

from PhysicsModules.AstroNeo.astro_neo.helper import AstroLogger, banner, bcolors
from PhysicsModules.AstroNeo.astro_neo.model import PARAM_NAMES
from PhysicsModules.AstroNeo.astro_neo.parser import parse_input_file
from PhysicsModules.AstroNeo.astro_neo.problem import AstroNeoProblem


class AstroNeo:
    """Top-level Astro Neo run: read -> setup -> run, mirroring
    ExafsNeo/NanoNeo's lifecycle. `astro_read` loads the .ini or a parameter
    dict, `astro_setup` builds the AstroNeoProblem and solver, `run` drives
    the generation loop and writes per-generation outputs."""

    def __init__(self, verbose_lvl: int = 5):
        self.logger = AstroLogger()
        self.verbose_lvl = verbose_lvl
        self.input_parameters = None
        self.problem = None
        self.solver = None
        self.output_path = None
        self.data_path = None
        self._gen_time = 0.0

    def astro_read(self, filepath: Union[str, pathlib.Path] = None, input_parameters: Dict = None) -> None:
        """Read the input file (or a parameter dict) into Astro Neo parameters."""
        if filepath is not None:
            file_path = os.path.join(os.getcwd(), str(filepath))
            self.input_parameters = parse_input_file(file_path)

        if input_parameters is not None:
            self.input_parameters = input_parameters

        if self.input_parameters is None:
            raise ValueError("No input parameters are given")

    def astro_setup(self) -> None:
        """Setup the problem, output files, and solver."""
        pars = self.input_parameters

        self._initialize_outputs(pars["output_file"])
        self.logger.print(banner())

        self.problem = AstroNeoProblem(
            data_dir=pars["data_dir"],
            data_file=pars["data_file"],
            bg_file=pars.get("bg_file"),
            rsp_file=pars.get("rsp_file"),
            acx2_path=pars.get("acx2_path"),
            xmin=pars.get("xmin", 7.0),
            xmax=pars.get("xmax", 30.0),
        )

        solver_cls = get_solver(pars.get("solOpt", 2))
        self.solver = solver_cls(
            self.problem,
            options={
                "nPops": pars.get("nPops", 20),
                "nGen": pars.get("nGen", 20),
                "selOpt": 0,
                "croOpt": pars.get("croOpt", 6),
                "mutOpt": pars.get("mutOpt", 4),
                "mutChance": pars.get("mutChance", 0.3),
                "nBestSample": pars.get("nBestSample", 0.3),
                "nLuckSample": pars.get("nLuckSample", 0.2),
                "F": pars.get("F", 0.5),
                "CR": pars.get("CR", 0.9),
            },
            logger=self.logger,
        )

    def run(self) -> SolverResult:
        """Run the fit; owns the generation loop for per-generation
        logging/output (same pattern as EXAFS/NanoIndentation)."""
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
        return solver.result

    # ------------------------- internals -------------------------

    def _initialize_outputs(self, output_file: str) -> None:
        self.output_path = os.path.join(os.getcwd(), output_file)
        self._check_if_exists(self.output_path)
        with open(self.output_path, "a+") as f:
            f.write("Gen,TPS,CURRFIT,BESTFIT\n")

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
        log(f"{bcolors.BOLD}Data File{bcolors.ENDC}: {self.problem.data_file}")
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
            log(f"{bcolors.BOLD}History Best: {bcolors.OKBLUE}{state.globBestVal}{bcolors.ENDC}")
        if hasattr(solver, "mutator"):
            log(f"Mutation Chance: {solver.mutator.mut_chance:.4f}")
        log("Time: " + str(round(self._gen_time, 5)) + "s")

    def _run_verbose_end(self, solver) -> None:
        log = self.logger.print
        log("-----------Output Stats---------------")
        log(f"{bcolors.BOLD}Total Time(s){bcolors.ENDC}: {round(solver.state.tt, 4)}")
        log(f"{bcolors.BOLD}Final Fitness{bcolors.ENDC}: {solver.state.globBestVal}")
        log("-------------------------------------------")

    def _output_generations(self, solver) -> None:
        """Output generation results into the run file and the best-fit data file."""
        state = solver.state
        with open(self.output_path, "a") as f1:
            f1.write(f"{state.currGen},{self._gen_time},{state.currBestVal},{state.globBestVal}\n")
        with open(self.data_path, "a") as f2:
            writer = csv.writer(f2)
            for name, value in zip(PARAM_NAMES, np.asarray(state.globBestInd.genes)):
                writer.writerow((name, value))
            f2.write("#################################\n")
