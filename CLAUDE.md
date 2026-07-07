# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository layout

NEO is a modular scientific computing framework:

- `Solvers/` — physics-agnostic optimization algorithms (numpy-only; **must never import larch or any physics module**). `core/` holds the framework (`ParameterSpace`/`GeneRange`, `Individual`, `Population`, `OptimizationProblem`, `RunState`, `SolverResult`, `BaseSolver`); `ga/` the genetic algorithm operators and `GASolver`/`GARechenbergSolver`; `de/` a Differential Evolution stub. `Solvers/__init__.py` has the solver registry (`get_solver("GA")`, numeric IDs 0/1/2 preserved from historical `solOpt` values).
- `PhysicsModules/EXAFS/` — fits EXAFS spectra with a GA (needs xraylarch).
- `PhysicsModules/NanoIndentation/` — Nano Neo: fits Oliver-Pharr power laws to nanoindentation unloading curves (numpy/matplotlib only). Ported from the standalone `nano_indent` package; its embedded GA was replaced by `Solvers`. Genome: `[A, hf, m] per path`, no shared genes. `nanoindentation_neo/gui/` is the legacy tkinter GUI, copied as-is and not wired to entry points.
- `PhysicsModules/XPS/` — XPS Neo: fits XPS spectra (Voigt/Gaussian/Lorentzian/Double Lorentzian/Doniach-Sunjic peaks, Baseline/Linear/Shirley/SVSC/Tougaard backgrounds) with its own GA/DE, ported from the standalone `XPS_Neo` package. **Does not use `Solvers`** — its genome is a heterogeneous list of floats interleaved with peak/background type-name strings, not a fixed-width float vector; see `PhysicsModules/XPS/README.md` ("Why this module doesn't use Solvers") before trying to route it through `Solvers.core`. Tests are **pytest**, not `unittest` — a deliberate exception, kept to preserve its bit-exact golden-master suite as-is.
- `PhysicsModules/common/` — non-physics infra shared across modules: `colors.py:TermColors` (ANSI codes) and `cli.py:str_to_bool`/`BaseLogger`. Each module's `helper.py` re-exports these under its historical names (`Bcolors`/`bcolors`, `NeoLogger`/`NanoLogger`); add new cross-module (non-solver) infra here rather than redefining it per module.

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
XPS tests (pytest, not unittest; needs numba/scipy):
```
cd PhysicsModules/XPS && pytest -m "not golden" -q     # unit + component, ~2s
cd PhysicsModules/XPS && pytest -m golden -q           # full CLI matrix, ~1min
```
Golden tests are bit-exact only in the pinned `PhysicsModules/XPS/constraints.txt` environment; in this repo's default env expect numpy-repr-only and last-digit float diffs (see `PhysicsModules/XPS/README.md`), not real regressions.

CLIs (entry points in `pyproject.toml`): `exafs_neo -i <input.ini>` (GUI via `exafs_neo_gui`), `nano_neo -i <input.ini>`, `xps_neo -i <input.ini>` (GUI via `xps_neo_gui`). Note the bundled EXAFS `tests/cu_test_files/test_cu.ini` references `path_files/Cu/...` paths that don't exist in this repo, so it is not runnable as-is; the actual test data lives at `PhysicsModules/EXAFS/tests/cu_test_files/cu_paths/`.

## Architecture

### Solvers ↔ physics coupling

The genome is a flat vector of genes; each gene samples from a discrete `GeneRange`. For EXAFS the layout is `[e0, (s02, sigma2, deltaR) per FEFF path]` — one shared E0 gene, three genes per path (defined in `exafs_neo/individual.py:build_parameter_space`). Generic operators only touch `pops.population`, `pops.next_population`, `pops.generate_individual()`, `individual.genes`, `pops.problem.fitness(genes)`, and `pops.state` (generation counter / best-fit tracking).

### EXAFS module flow

`input_arg.py` (CLI) → `parser.py`/`ini_parser.py` (.ini → validated dict) → `neoPars.py:NeoPars` (central config/state, attrs-based sub-par classes; `.ini` keys map here) → `exafs.py:ExafsNeo` runs the loop: per generation `NeoSolver.solve(...)` → selection/crossover/mutation → `NeoPopulations.eval_population()` → `NeoResult.collect` + `NeoFilePars` output writing. E0 gets specially optimized at `nGen//2` and at the end (`exafs_pop.py:optimize_e0`).

**Back-compat shim layer**: `neoSelector.py`, `neoCrossOver.py`, `neoMutator.py`, `neoSolver.py` keep the historical EXAFS API (`NeoSelector().initialize(exafs_pars)` / `.select(pops)` etc.) but delegate the algorithms to `Solvers.ga`. Preserved quirks the tests pin: `mutator.mutOpt == mut_options + 1`, exact `mutType`/`croType`/`solver_operator` label strings, and Rechenberg writing `mutPars.mutChance` (which the active mutator sampled at init — historically a near-no-op). Don't "fix" these without updating the tests.

`exafs_neo/problem.py` also provides `NeoRunStateView`, a read-only adapter exposing the Solvers `RunState` interface over `NeoPars` bookkeeping (`runPars`/`bestFitPars`), attached to `NeoPopulations` as `pops.state` alongside `pops.problem`.

### XPS module flow

`input_arg.py` (CLI, argparse + `--seed`) → `parser2.py`/`ini_parser.py` (`.ini` → flat `config` dict via `load_config`) → `xps.py:main()` does `globals().update(config)` — `XPS_GA`'s methods deliberately read config as bare module globals of `xps.py` (a documented transitional seam carried over from the source project; don't "fix" a bare read into a `self.x`/`cfg.x` write, that changes mutation-sharing semantics for `addPeak`/`removePeak`). `XPS_GA.run()` drives the GA/DE loop directly (its own `next_generation`/`mutatePopulation`/`crossover`, or the DE path when `mutated_options == 3`) — no `Solvers` involvement; see `PhysicsModules/XPS/README.md` for why. `xps_individual.py:Individual.getFit` evaluates a candidate by summing `xps_fit.py` peak/background curves; `loss.py:compute_loss` is the one fitness function, shared with the GUI's post-analysis.

### Adding a new solver

Implement `step()` on a `BaseSolver` subclass under `Solvers/<name>/`, register it in `SOLVER_REGISTRY` in `Solvers/__init__.py`. The base class owns the generation loop, result collection, and problem hooks.

### Adding a new physics module

Two patterns, pick based on whether the fit parameters are a fixed-width vector of floats:
- **Flat float genome** (most cases): subclass `Solvers.core.OptimizationProblem` (all I/O inside the module), add an input parser + CLI entry point in `pyproject.toml` `[project.scripts]`, and module tests under `PhysicsModules/<Name>/tests/` (`unittest`, run via `python -m unittest discover`). `PhysicsModules/EXAFS/exafs_neo/problem.py` and `PhysicsModules/NanoIndentation/nanoindentation_neo/problem.py` are the reference implementations.
- **Heterogeneous/variable-shape genome**: keep the module's own solver loop self-contained, as `PhysicsModules/XPS/` does — don't force-fit it through `Solvers.core`. Document the decision in the module's README the way `PhysicsModules/XPS/README.md` does.

Either way: add the test run to both `.github/workflows/*.yml`.

## Gotchas

- `attrs` classes in `neoPars.py`: use `field(factory=list)`, never `field(default=[])` — a shared mutable default here previously made `pathrange_pars` grow across `NeoPars` instances in one process.
- EXAFS mutation chances are historically **percent-scaled** in some paths (`random() * 100 < chance` in `Individual.mutate_paths` and `mutate_e0`) but probability-scaled in others (`random() < mutChance` for per-individual selection). The EXAFS `Individual.mutate()` override deliberately keeps percent semantics; generic `Solvers.core.Individual.mutate()` is probability-based.
- Self-absorption correction is an optional git submodule (`contrib/sabcor`) built with `make`; not part of the default install.
