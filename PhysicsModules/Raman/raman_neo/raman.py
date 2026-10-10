"""
Raman Neo runner: wires the RamanProblem to the shared NEO solvers.

Mirrors the EXAFS/NanoIndentation module pattern (read -> setup -> run); the
GA/DE itself lives in Solvers. Historical .ini option IDs map straight onto
Solvers:

    solver_options    -> get_solver(...)  0 GA, 1 GA + Rechenberg, 2 DE
                         (or any registry name: PSO, CMA_ES, ...)
    mutated_options   -> Solvers mutOpt   0 per individual, 1 per gene
                         ("per band"), 2 per trait (stub), 3 Metropolis,
                         4 bounded
    selection_options -> Solvers selOpt;  crossover_options -> Solvers croOpt
"""

import csv
import os
import pathlib
import time
from typing import Dict, Union

import numpy as np

from Solvers import get_solver
from Solvers.core import SolverResult

from PhysicsModules.Raman.raman_neo.analysis.raman_analysis import RamanAnalysis
from PhysicsModules.Raman.raman_neo.helper import RamanLogger, banner, bcolors
from PhysicsModules.Raman.raman_neo.parser import parse_input_file
from PhysicsModules.Raman.raman_neo.problem import RamanProblem


class RamanNeo:
    """Top-level Raman Neo run: read -> setup -> run, mirroring
    ExafsNeo/NanoNeo's lifecycle. `raman_read` loads the .ini or a parameter
    dict, `raman_setup` builds the RamanProblem and solver, `run` drives the
    generation loop and writes per-generation outputs."""

    def __init__(self, verbose_lvl: int = 5):
        self.logger = RamanLogger()
        self.verbose_lvl = verbose_lvl
        self.input_parameters = None
        self.problem = None
        self.solver = None
        self.output_path = None
        self.data_path = None
        self.peaks_path = None
        self.band_optimize_result = None
        self.analysis = None
        self._gen_time = 0.0

    def raman_read(self, filepath: Union[str, pathlib.Path] = None, input_parameters: Dict = None) -> None:
        """Read the input file (or a parameter dict) into Raman Neo parameters."""
        if filepath is not None:
            file_path = os.path.join(os.getcwd(), str(filepath))
            self.input_parameters = parse_input_file(file_path)

        if input_parameters is not None:
            self.input_parameters = input_parameters

        if self.input_parameters is None:
            raise ValueError("No input parameters are given")

    def raman_setup(self) -> None:
        """Setup the problem, output files, and solver."""
        pars = self.input_parameters

        self._initialize_outputs(pars["output_file"])
        if self.verbose_lvl > 0:
            self.logger.print(banner())

        self.problem = RamanProblem(
            data_file=pars["data_file"],
            band_file=pars["band_file"],
            band_list=pars.get("band_list"),
            xmin=pars.get("xmin", 100.0),
            xmax=pars.get("xmax", 1800.0),
            deltax=pars.get("deltax", 1.0),
            weight=pars.get("weight", 0.0),
            baseline_order=pars.get("baseline_order", 0),
            normalize=pars.get("normalize", True),
            bandrange_file=pars.get("bandrange_file"),
            num_compounds=pars.get("num_compounds", 1),
            shift_range=pars.get("shift_range", ()),
            shift_sweep=pars.get("shiftSweep", True),
        )

        solver_cls = get_solver(pars.get("solOpt", 0))
        self.solver = solver_cls(
            self.problem,
            options={
                "nPops": pars.get("nPops", 100),
                "nGen": pars.get("nGen", 100),
                "selOpt": pars.get("selOpt", 0),
                "croOpt": pars.get("croOpt", 0),
                "mutOpt": pars.get("mutOpt", 1),
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
        logging/output (same order as Solvers' BaseSolver.run, including the
        problem's on_generation_end/on_run_end shift-sweep hooks)."""
        solver = self.solver
        problem = self.problem
        self._run_verbose_start()

        solver.initialize()
        for _ in range(solver.n_gen):
            solver.state.start_gen()
            gen_start = time.time()
            solver.step()
            solver.result.collect(solver.state, solver.population)
            problem.on_generation_end(solver.state, solver.population)
            self._gen_time = time.time() - gen_start
            self._run_verbose_gen(solver)
            self._output_generations(solver)
            solver.state.end_gen()
        problem.on_run_end(solver.state, solver.population)
        # the final shift sweep can only improve the (same) best individual
        solver.result.best_individual = solver.state.globBestInd
        solver.result.best_value = solver.state.globBestVal

        best_genes = solver.state.globBestInd.genes
        if self.input_parameters.get("bandOptimize", False):
            self.band_optimize_result = problem.band_optimize(
                best_genes, self.input_parameters.get("bandOptimizePercent", 0.01)
            )
            self._run_verbose_band_optimize()

        self.analysis = RamanAnalysis(problem, best_genes, fitness=solver.state.globBestVal)
        self._write_peak_table()
        self._run_verbose_end(solver)
        if self.input_parameters.get("printGraph", False):
            self.plot(solver.state.globBestInd)
        return solver.result

    def plot(self, individual) -> None:
        import matplotlib.pyplot as plt

        RamanAnalysis(self.problem, individual.genes).plot(show_bands=True)
        plt.show()

    # ------------------------- internals -------------------------

    def _initialize_outputs(self, output_file: str) -> None:
        self.output_path = os.path.join(os.getcwd(), output_file)
        self._check_if_exists(self.output_path)
        with open(self.output_path, "a+") as f:
            f.write("Gen,TPS,FITTNESS,CURRFIT,CURRIND,BESTFIT,BESTIND\n")

        base = os.path.splitext(self.output_path)[0]
        self.data_path = base + "_data.csv"
        self._check_if_exists(self.data_path)
        self.peaks_path = base + "_peaks.csv"
        self._check_if_exists(self.peaks_path)

        log_path = base + ".log"
        self._check_if_exists(log_path)
        self.logger.initialize_logging(log_path)

    @staticmethod
    def _check_if_exists(path_file) -> None:
        if os.path.exists(path_file):
            os.remove(path_file)
        pathlib.Path(path_file).parent.mkdir(parents=True, exist_ok=True)

    def _run_verbose_start(self) -> None:
        if self.verbose_lvl <= 0:
            return
        log = self.logger.print
        problem = self.problem
        log("-----------Inputs File Stats---------------")
        log(f"{bcolors.BOLD}File{bcolors.ENDC}: {problem.data_file}")
        log(f"{bcolors.BOLD}Band File{bcolors.ENDC}: {problem.band_file}")
        log(f"{bcolors.BOLD}Bands{bcolors.ENDC}: {problem.band_names}")
        log(f"{bcolors.BOLD}Window (cm^-1){bcolors.ENDC}: [{problem.spectrum.xmin}, {problem.spectrum.xmax}]")
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
        with np.printoptions(precision=4, suppress=True):
            log(f"{bcolors.BOLD}Best fit: {bcolors.OKBLUE}{state.currBestVal}{bcolors.ENDC}")
            log("Best fit combination [amp, shift, width, center_shift]:\n"
                + str(self.problem.band_params(state.currBestInd.genes)))
            log(f"{bcolors.BOLD}History Best: {bcolors.OKBLUE}{state.globBestVal}{bcolors.ENDC}")
        mutator = getattr(solver, "mutator", None)
        if mutator is not None:
            log(f"Mutation Chance: {mutator.mut_chance:.4f}")
        log("Time: " + str(round(self._gen_time, 5)) + "s")

    def _run_verbose_band_optimize(self) -> None:
        if self.verbose_lvl <= 0:
            return
        log = self.logger.print
        result = self.band_optimize_result
        log("-----------Band Optimize---------------")
        log(f"Threshold: {100 * result['band_optimize_percent']:.2f}% of total band area")
        for name, fraction in result["optimize_band_lists"].items():
            log(f"  {name}: {100 * fraction:.2f}%")
        dropped = [b for b in result["original_bands"] if b not in result["optimize_bands"]]
        if dropped:
            log(f"Below threshold: {dropped}")

    def _run_verbose_end(self, solver) -> None:
        if self.verbose_lvl <= 0:
            return
        log = self.logger.print
        log("-----------Output Stats---------------")
        log(f"{bcolors.BOLD}Total Time(s){bcolors.ENDC}: {round(solver.state.tt, 4)}")
        log(f"{bcolors.BOLD}Final Fitness{bcolors.ENDC}: {solver.state.globBestVal}")
        log(self.analysis.report())

    def _output_generations(self, solver) -> None:
        """Output generation results into the run file and the best-fit data
        file ([amp, shift, width, center_shift] per band, EXAFS layout)."""
        state = solver.state
        with open(self.output_path, "a") as f1:
            f1.write(
                f"{state.currGen},{self._gen_time},"
                f"{state.currBestVal},{list(state.currBestInd.genes)},"
                f"{state.globBestVal},{list(state.globBestInd.genes)}\n"
            )
        with open(self.data_path, "a") as f2:
            writer = csv.writer(f2)
            for row in self.problem.band_params(state.globBestInd.genes):
                writer.writerow(tuple(row))
            f2.write("#################################\n")

    def _write_peak_table(self) -> None:
        self.analysis.write_peak_table(self.peaks_path)
