# NEO Framework Design

## Motivation

NEO began as EXAFS Neo: a genetic algorithm hard-wired to EXAFS spectrum
fitting, with the GA operators, population bookkeeping, and physics all in
one package. The same GA was copy-pasted into sibling projects (nano_indent,
and prospectively XPS), each fork drifting independently.

The goal of this design is a modular framework where:

- **optimization algorithms live in exactly one place** (`Solvers/`) and are
  physics-agnostic;
- **each physics application is a self-contained module** under
  `PhysicsModules/`, owning its data formats, model, input parsing, CLI, and
  tests;
- a new physics module plugs into every existing solver — and a new solver
  becomes available to every physics module — **without changes on the other
  side**.

## Architecture

```
┌──────────────────────────────────────────────────────────┐
│ PhysicsModules/                                          │
│                                                          │
│  EXAFS/          NanoIndentation/       XPS/             │
│  (larch, FEFF)   (Oliver-Pharr)         (scaffold)       │
│      │                 │                  │              │
│      └────────┬────────┴──────────────────┘              │
│               ▼                                          │
│    OptimizationProblem (the contract)                    │
│    - space: ParameterSpace                               │
│    - fitness(genes) -> float        (lower = better)     │
│    - sample_genes()                 (optional bias)      │
│    - on_generation_end / on_run_end (optional hooks)     │
└───────────────┬──────────────────────────────────────────┘
                ▼
┌──────────────────────────────────────────────────────────┐
│ Solvers/  (numpy only — never imports physics)           │
│                                                          │
│  core/  ParameterSpace · GeneRange · Individual          │
│         Population · RunState · SolverResult             │
│         BaseSolver (owns the generation loop)            │
│  ga/    selectors · crossovers · mutators                │
│         GASolver · GARechenbergSolver                    │
│  de/    DESolver (stub)                                  │
│                                                          │
│  registry: get_solver("GA" | "GA_Rechenberg" | "DE")     │
└──────────────────────────────────────────────────────────┘
```

### The genome model

A candidate solution is a **flat vector of genes**. Each gene draws from a
`GeneRange`: a discrete, ordered grid of allowed values — matching how the
physics packages have always defined fitting-parameter grids (`np.arange`
over low/high/step). A `ParameterSpace` is the ordered list of gene ranges
and defines sampling, clipping, and per-gene limits.

Physics modules define the genome layout and its meaning:

| Module          | Genome layout                              | Shared genes |
|-----------------|--------------------------------------------|--------------|
| EXAFS           | `[e0, (s02, sigma2, deltaR) × npath]`      | e0           |
| NanoIndentation | `[(A, hf, m) × npaths]`                    | none         |
| XPS (future)    | e.g. `[(BE, amp, fwhm, mix) × npeaks]`     | tbd          |

Generic operators never interpret genes — they only sample, mix, and perturb
them within the space. Anything meaning-dependent (EXAFS's shared-E0
mutation, its mid-run E0 grid sweep) belongs to the physics module, expressed
through the problem hooks or module-level orchestration.

### The solver side

`BaseSolver` owns the run: initialize population → per generation call
`step()` (implemented by each solver), collect into `SolverResult`, fire
`problem.on_generation_end`, and finally `problem.on_run_end`. Operators are
pluggable strategies selected by the same numeric option IDs the historical
`.ini` files used (`selOpt`, `croOpt`, `mutOpt`, `solOpt`):

- **Selectors**: roulette wheel (0), tournament (1, stub)
- **Crossovers**: uniform (0), single point (1), dual point (2, stub),
  arithmetic (3), or (4), average (5)
- **Mutators**: per individual (0), per gene (1), per trait (2, stub),
  Metropolis (3), bounded (4)
- **Solvers**: GA (0), GA + Rechenberg 1/5-rule mutation adaptation (1),
  DE (2, stub)

`RunState` carries generation counters, timing, and best-fit tracking, so
operators that need run context (Metropolis temperature, Rechenberg counter,
bounded-step shrinkage) read it from one place instead of reaching into
physics config objects.

### Module anatomy (the contract)

Every module under `PhysicsModules/<Name>/` provides:

1. `problem.py` — the `OptimizationProblem` subclass. All I/O and model math
   stay behind this interface; solvers never see files or spectra.
2. An input parser for the module's historical `.ini` format (functions,
   not import-time side effects) and a CLI entry point registered in the
   root `pyproject.toml` under `[project.scripts]`.
3. `tests/` — run from the repository root; must not require instrument
   data when synthetic fixtures suffice.
4. `README.md` — setup, input format, examples.

`PhysicsModules/XPS/README.md` documents this contract for the next module.

## Module status and porting notes

### EXAFS (functional; back-compat frozen)

The original public API was kept intact: same `exafs_neo` CLI, same `.ini`
keys, same test suite (verified identical pass/fail against the pre-refactor
baseline). The operator modules (`neoSelector`, `neoCrossOver`, `neoMutator`,
`neoSolver`) are now thin shims delegating to `Solvers.ga`, preserving pinned
quirks (e.g. `mutOpt == mut_options + 1`, exact label strings, Rechenberg
writing `mutPars.mutChance`). `NeoRunStateView` adapts `NeoPars` bookkeeping
to the `RunState` interface so generic operators see live EXAFS state.
EXAFS's `Individual` subclasses the Solvers individual over the flat genome
while keeping the historical path-oriented accessors.

### NanoIndentation (functional; clean re-port)

Ported from the standalone `nano_indent` package with no API freeze
(nothing depended on its internals). Its embedded GA — a fork of the EXAFS
one — was deleted in favor of `Solvers`; the module keeps only physics:
`NanoIndent_Data` (Hysitron/iMicro/CSV loading + unloading-segment slicing),
the Oliver-Pharr model, `.ini` parsing (same file format), and the `NanoNeo`
runner (`nano_neo` CLI) with the historical per-generation output files.
Deliberate behavior changes: NaN fitness → `inf` (instead of
discard-and-regenerate), import-time argparse/globals removed, Rechenberg
step unified with the shared schedule. The old tkinter GUI is carried as-is
under `nanoindentation_neo/gui/`, unwired.

### XPS (scaffold)

Placeholder `XPSProblem` (fitness raises `NotImplementedError`), contract
tests proving it plugs into the registry, README describing the fill-in
steps.

## Design rules

1. `Solvers/` imports numpy only. A physics dependency appearing there is a
   regression — it breaks the fast, dependency-free solver test suite.
2. Physics modules never import each other.
3. Historical numeric option IDs and `.ini` keys are stable interfaces;
   changing them requires migrating the affected module's input format
   deliberately, not as a side effect.
4. New solver = subclass `BaseSolver`, implement `step()`, register in
   `SOLVER_REGISTRY`. New physics = implement `OptimizationProblem`, add
   parser + CLI + tests + README, add the test path to CI.

## Verification strategy

- `Solvers/tests/` proves the framework on an analytic toy problem
  (convergence on `sum((x - target)^2)`), with no physics installed.
- EXAFS parity was established by running the full suite from the
  pre-refactor commit and after the refactor in the same larch environment:
  identical results (89 tests; 5 pre-existing `test_EXAFS_analysis` errors),
  plus an end-to-end GA run on the bundled Cu data producing comparable
  fitness trajectories.
- NanoIndentation is verified against synthetic Oliver-Pharr curves with
  known ground truth (recovered within grid tolerance), plus an end-to-end
  `.ini` → outputs run for every mutation option.
