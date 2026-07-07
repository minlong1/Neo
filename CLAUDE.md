# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository layout

NEO is a modular scientific computing framework:

- `Solvers/` — physics-agnostic optimization algorithms (numpy-only; **must never import larch or any physics module**). `core/` holds the framework (`ParameterSpace`/`GeneRange`, `Individual`, `Population`, `OptimizationProblem`, `RunState`, `SolverResult`, `BaseSolver`); `ga/` the genetic algorithm operators and `GASolver`/`GARechenbergSolver`; `de/` a Differential Evolution stub. `Solvers/__init__.py` has the solver registry (`get_solver("GA")`, numeric IDs 0/1/2 preserved from historical `solOpt` values).
- `PhysicsModules/EXAFS/` — fits EXAFS spectra with a GA (needs xraylarch).
- `PhysicsModules/NanoIndentation/` — Nano Neo: fits Oliver-Pharr power laws to nanoindentation unloading curves (numpy/matplotlib only). Ported from the standalone `nano_indent` package; its embedded GA was replaced by `Solvers`. Genome: `[A, hf, m] per path`, no shared genes. `nanoindentation_neo/gui/` is the legacy tkinter GUI, copied as-is and not wired to entry points.
- `PhysicsModules/XPS/` — scaffold only; placeholder `problem.py` with `fitness()` raising `NotImplementedError`, contract tests, and a README describing what to fill in.

A physics module plugs into the solvers by subclassing `Solvers.core.OptimizationProblem` (a `ParameterSpace` + `fitness(genes) -> float`, plus optional `sample_genes`/`on_generation_end`/`on_run_end` hooks). `PhysicsModules/EXAFS/exafs_neo/problem.py:EXAFSProblem` is the reference implementation.

## Commands

Run all commands from the **repository root** unless noted — the codebase uses absolute imports rooted there (e.g. `from PhysicsModules.EXAFS.exafs_neo.exafs import ExafsNeo`).

Install (editable):
```
pip install -e .
```
`xraylarch` is finicky — see `PhysicsModules/EXAFS/README.md` (conda/mamba recommended). The `Solvers` package and its tests need only numpy.

Solvers tests (fast, no larch needed):
```
python -m unittest discover -s Solvers/tests -t . -v
```

EXAFS tests (needs xraylarch; **cwd must be the EXAFS module** because test data paths like `tests/cu_test_files/...` are cwd-relative, while `-t ../..` roots imports at the repo):
```
cd PhysicsModules/EXAFS && python -m unittest discover -s tests -t ../.. -v
```
Equivalent: `PhysicsModules/EXAFS/run_tests`. Single test:
```
cd PhysicsModules/EXAFS && python -m unittest PhysicsModules.EXAFS.tests.test_neosolver -v   # needs repo root on sys.path, e.g. PYTHONPATH=../..
```
Known baseline: the 5 tests in `test_EXAFS_analysis.py` error (broken data paths, pre-existing). Everything else passes.

NanoIndentation tests (numpy-only; synthetic data, no instrument files needed):
```
python -m unittest discover -s PhysicsModules/NanoIndentation/tests -t . -v
```
XPS scaffold tests: `python -m unittest discover -s PhysicsModules/XPS/tests -t . -v`.

CLIs (entry points in `pyproject.toml`): `exafs_neo -i <input.ini>` (GUI via `exafs_neo_gui`), `nano_neo -i <input.ini>`. Note the bundled EXAFS `tests/cu_test_files/test_cu.ini` references `path_files/Cu/...` paths that don't exist in this repo, so it is not runnable as-is; the actual test data lives at `PhysicsModules/EXAFS/tests/cu_test_files/cu_paths/`.

## Architecture

### Solvers ↔ physics coupling

The genome is a flat vector of genes; each gene samples from a discrete `GeneRange`. For EXAFS the layout is `[e0, (s02, sigma2, deltaR) per FEFF path]` — one shared E0 gene, three genes per path (defined in `exafs_neo/individual.py:build_parameter_space`). Generic operators only touch `pops.population`, `pops.next_population`, `pops.generate_individual()`, `individual.genes`, `pops.problem.fitness(genes)`, and `pops.state` (generation counter / best-fit tracking).

### EXAFS module flow

`input_arg.py` (CLI) → `parser.py`/`ini_parser.py` (.ini → validated dict) → `neoPars.py:NeoPars` (central config/state, attrs-based sub-par classes; `.ini` keys map here) → `exafs.py:ExafsNeo` runs the loop: per generation `NeoSolver.solve(...)` → selection/crossover/mutation → `NeoPopulations.eval_population()` → `NeoResult.collect` + `NeoFilePars` output writing. E0 gets specially optimized at `nGen//2` and at the end (`exafs_pop.py:optimize_e0`).

**Back-compat shim layer**: `neoSelector.py`, `neoCrossOver.py`, `neoMutator.py`, `neoSolver.py` keep the historical EXAFS API (`NeoSelector().initialize(exafs_pars)` / `.select(pops)` etc.) but delegate the algorithms to `Solvers.ga`. Preserved quirks the tests pin: `mutator.mutOpt == mut_options + 1`, exact `mutType`/`croType`/`solver_operator` label strings, and Rechenberg writing `mutPars.mutChance` (which the active mutator sampled at init — historically a near-no-op). Don't "fix" these without updating the tests.

`exafs_neo/problem.py` also provides `NeoRunStateView`, a read-only adapter exposing the Solvers `RunState` interface over `NeoPars` bookkeeping (`runPars`/`bestFitPars`), attached to `NeoPopulations` as `pops.state` alongside `pops.problem`.

### Adding a new solver

Implement `step()` on a `BaseSolver` subclass under `Solvers/<name>/`, register it in `SOLVER_REGISTRY` in `Solvers/__init__.py`. The base class owns the generation loop, result collection, and problem hooks.

### Adding a new physics module

Follow `PhysicsModules/XPS/README.md`: subclass `OptimizationProblem` (all I/O inside the module), add an input parser + CLI entry point in `pyproject.toml` `[project.scripts]`, and module tests under `PhysicsModules/<Name>/tests/`. Add the test run to both `.github/workflows/*.yml`.

## Gotchas

- `attrs` classes in `neoPars.py`: use `field(factory=list)`, never `field(default=[])` — a shared mutable default here previously made `pathrange_pars` grow across `NeoPars` instances in one process.
- EXAFS mutation chances are historically **percent-scaled** in some paths (`random() * 100 < chance` in `Individual.mutate_paths` and `mutate_e0`) but probability-scaled in others (`random() < mutChance` for per-individual selection). The EXAFS `Individual.mutate()` override deliberately keeps percent semantics; generic `Solvers.core.Individual.mutate()` is probability-based.
- Self-absorption correction is an optional git submodule (`contrib/sabcor`) built with `make`; not part of the default install.
