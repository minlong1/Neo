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
│  EXAFS/          NanoIndentation/       XPS/*            │
│  (larch, FEFF)   (Oliver-Pharr)      (own GA/DE loop)    │
│      │                 │                                 │
│      └────────┬────────┘                                 │
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
│  de/    DESolver · differential_evolution_step()         │
│  demcmc/ DEMCMCSolver — DE-BNN (Bayesian NN via DE-MCMC) │
│                                                          │
│  registry: get_solver("GA"|"GA_Rechenberg"|"DE"|"DE_MCMC")│
└──────────────────────────────────────────────────────────┘
```

\* XPS doesn't route through `OptimizationProblem`/`Solvers` — its genome
is a heterogeneous, type-tagged list rather than a fixed-width float
vector, so it keeps its own (already working, already tested) GA/DE loop.
See "Module status and porting notes" below and
`PhysicsModules/XPS/README.md`.

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
| XPS             | not a flat float vector — see below        | n/a          |

XPS's actual "genome" (`Individual.get_params()`) is a list of floats
**interleaved with peak/background type-name strings** (`'Voigt'`,
`'Baseline'`, ...) marking block boundaries, and the number of floats per
block varies with peak type (singlet vs. doublet vs. Coster-Kronig) and
correlated/limited-parameter settings. That doesn't fit the
`ParameterSpace` model above, so XPS opts out of `Solvers` entirely rather
than distorting either side to make it fit.

Generic operators never interpret genes — they only sample, mix, and perturb
them within the space. Anything meaning-dependent (EXAFS's shared-E0
mutation, its mid-run E0 grid sweep) belongs to the physics module, expressed
through the problem hooks or module-level orchestration.

Not every genome is a discrete grid, though: `Solvers.demcmc` (DE-BNN, see
below) trains neural-network weights, which are continuous and unbounded
(He/Xavier-style initialization). `ContinuousGeneRange` is the sibling to
`GeneRange` for this case — same interface (`sample`/`clip`/`low`/`high`),
sampled from `N(mean, std^2)` instead of a materialized value array — so
`ParameterSpace` and every generic operator work unchanged.

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
  DE/rand/1/bin (2)

`RunState` carries generation counters, timing, and best-fit tracking, so
operators that need run context (Metropolis temperature, Rechenberg counter,
bounded-step shrinkage) read it from one place instead of reaching into
physics config objects.

DE is also available decoupled from `BaseSolver`'s own generation loop, as
`Solvers.de.differential_evolution_step(population, F, CR)` — a physics
module that owns its own top-level loop (EXAFS's `ExafsNeo.run()`) calls
this directly, per generation, against its own population object, rather
than duplicating the DE math or handing its loop over to `DESolver.run()`.
The function only needs `population.problem` (`.space.clip`, `.fitness`),
`population.population` (a list of individuals exposing `.genes`),
`population.eval_population(replace, sorting)`, and
`population.generate_individual()` — both `Solvers.core.Population` and
EXAFS's `NeoPopulations` satisfy this without adapters.

### DE-BNN: differential evolution as MCMC (`Solvers/demcmc/`)

Forbes & Long, "DE-BNN: An evolutionary approach to Bayesian neural
network posterior sampling" (*Neurocomputing* 678 (2026) 133103), observe
that DE's greedy selection (accept a trial if it scores better) is already
an MCMC acceptance rule; adding small Gaussian noise to the mutant for
detailed balance turns the sequence of accepted candidates at each
population index into a Markov chain — DE gets `nPops` chains essentially
for free. Applied to training a neural network, the population of
post-burn-in samples *is* the posterior over weights/biases, enabling
Bayesian prediction (mean + credible interval) instead of a single trained
network. The paper's own conclusion names its intended home explicitly:
"future work… incorporate DE-MCMC and DE-BNN as an available solver in the
NEO platform" — i.e. this repository — which is why it lives in `Solvers/`
as a new solver rather than a `PhysicsModules` domain application; DE-BNN
is a general regression/training capability, not tied to a physics domain.

Implemented, mapped to the paper's sections:
- Base loop + configurable mutation operator (`rand/1`, `rand/2`, `best/1`,
  `best/2`, Eq. 5-8) + MCMC noise term (Eq. 48) — `mutation.py`,
  `demcmc_solver.py:de_mcmc_generation_step`.
- Hyper-mutation (Section 3.1-3.4): a running-average-residual stagnation
  detector resamples F/CR/operator (Eq. 11, 12, 15) while the search has
  stopped improving — `hyper_mutation.py`.
- SVD and local-search refinement (Section 3.5, 3.7, Eq. 16-25): run every
  `refine_every` generations, keep a modified candidate only if it beats
  the *population's* current best (the paper's stated criterion for both —
  a deliberately strict bar, not "improves over its own previous value") —
  `refinement.py`.
- Single- and multi-chain posterior collection and prediction (Section 5):
  `posterior.py:PosteriorResult` — `.predict()` averages forward passes
  over posterior samples (Eq. 49-50: the statistically correct predictive
  posterior, not a forward pass through averaged weights).

Not implemented: clustering-based refinement (Section 3.6, k-means/
spectral/agglomerative) — one of three interchangeable refinement
techniques, deliberately skipped rather than pulling in a clustering
dependency (e.g. scikit-learn) `Solvers` otherwise doesn't need.
`cluster_refine` is a documented `NotImplementedError` stub, matching this
repo's existing pattern for other intentionally-unimplemented variants
(`TournamentSelector`, `DualPointCrossover`, `PerTraitMutator` in
`Solvers/ga`).

Verified against the paper's own worked example (a (7,35,10,5,1) MLP has
701 parameters, Eq. 51-52) and end-to-end on synthetic sine regression:
loss decreases, posterior samples collect at the expected count
(`(nGen - burn_in) * num_chains`), and the predictive posterior mean
tracks the true function with the credible interval covering it.

### Module anatomy (the contract)

Every module under `PhysicsModules/<Name>/` provides:

1. Either `problem.py` — the `OptimizationProblem` subclass (all I/O and
   model math stay behind this interface; solvers never see files or
   spectra) — **or**, if the fit parameters genuinely aren't a fixed-width
   float vector (see XPS), a self-contained solver loop, with the decision
   and why documented in the module's README.
2. An input parser for the module's historical `.ini` format (functions,
   not import-time side effects) and a CLI entry point registered in the
   root `pyproject.toml` under `[project.scripts]`.
3. `tests/` — must not require instrument data when synthetic fixtures
   suffice. `unittest`, run from the repository root, for the `Solvers`-based
   modules; XPS is the one exception (pytest, kept as ported — see its
   README).
4. `README.md` — setup, input format, examples.

## Module status and porting notes

### EXAFS (functional; back-compat frozen)

The original public API was kept intact: same `exafs_neo` CLI, same `.ini`
keys, same test suite (verified identical pass/fail against the pre-refactor
baseline). The operator modules (`neoSelector`, `neoCrossOver`, `neoMutator`,
`neoSolver`) are now thin shims delegating to `Solvers.ga`/`Solvers.de`,
preserving pinned quirks (e.g. `mutOpt == mut_options + 1`, exact label
strings, Rechenberg writing `mutPars.mutChance`). `NeoRunStateView` adapts
`NeoPars` bookkeeping to the `RunState` interface so generic operators see
live EXAFS state. EXAFS's `Individual` subclasses the Solvers individual
over the flat genome while keeping the historical path-oriented accessors.
`solOpt == 2` (`NeoSolver_DE`, previously an unimplemented stub matching
upstream EXAFS_Neo) now delegates to `Solvers.de.differential_evolution_step`
directly against `NeoPopulations` — the first of the three EXAFS solver
options to run actual differential evolution rather than a no-op.

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

### XPS (functional; ported as-is, self-contained GA)

Ported from the standalone `XPS_Neo` package, which had already gone
through its own production-readiness effort (determinism/seeding, a
golden-master regression harness, killing import-time side effects,
unifying a diverged GUI/CLI fork, a performance pass — all bit-exact
gated). That work is carried over essentially mechanically (import paths
rewritten to `PhysicsModules.XPS.xps_neo.*`; `helper.py` reuses
`PhysicsModules/common` for colors/`str_to_bool` the same way EXAFS and
NanoIndentation do). Unlike the other two modules, XPS's GA/DE loop
(`xps.py:XPS_GA`) is **not** rewired onto `Solvers` — see "The genome
model" above and `PhysicsModules/XPS/README.md` for the reasoning. Its
test suite (pytest: unit, component/numerical, and golden-master) is
likewise ported as-is rather than translated to `unittest`, to keep its
bit-exact regression discipline intact.

## Design rules

1. `Solvers/` imports numpy only. A physics dependency appearing there is a
   regression — it breaks the fast, dependency-free solver test suite.
2. Physics modules never import each other.
3. Historical numeric option IDs and `.ini` keys are stable interfaces;
   changing them requires migrating the affected module's input format
   deliberately, not as a side effect.
4. New solver = subclass `BaseSolver`, implement `step()`, register in
   `SOLVER_REGISTRY`. New physics = implement `OptimizationProblem` (or, if
   the genome genuinely isn't a fixed-width float vector, a self-contained
   solver loop, documented like XPS's) — add parser + CLI + tests + README,
   add the test path to CI.

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
- XPS parity was checked against its own golden-master matrix (9 cases
  spanning every peak/background family, singlet/doublet/Coster-Kronig,
  correlated/limited parameters, and both the GA and DE paths): all 9
  reproduce their expected generation-by-generation trajectory once the
  environment-dependent `np.float64(...)` repr difference (numpy ≥2.0 vs.
  this repo's numpy 1.26) is normalized away; 4 of 9 are bit-exact even
  without normalizing. The remaining divergence is confirmed last-significant-digit
  float noise from running outside the pinned `constraints.txt` environment
  (numpy/numba/scipy version drift) — not a logic difference. 34 of 36 fast
  (non-golden) tests pass outright; the other 2 fail at the ~1e-14 relative
  level for the same reason.
- DE-BNN (`Solvers/demcmc/`) is verified against the paper's own 701-
  parameter worked example (dimension count only, Eq. 51-52) and 26 unit/
  integration tests: MLP flatten/unflatten round-trips, all four mutation
  operators run, the stagnation tracker correctly flags a flat vs.
  improving fitness sequence, SVD/local-search refinement never worsen the
  population's best score, and a full `DEMCMCSolver` run on synthetic sine
  regression converges with posterior sample counts matching
  `(nGen - burn_in) * num_chains` and a credible interval that brackets
  the true function.
