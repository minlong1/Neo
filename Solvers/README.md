# Solvers

Physics-agnostic, population-based optimization algorithms for NEO. Pure
numpy — this package must never import `larch`, `xspec`, or any other
physics-module dependency, so it stays usable (and testable) without any of
them installed.

## The plug-in contract

A physics module exposes its fitting task as an `OptimizationProblem`: a
`ParameterSpace` describing the genome (one `GeneRange` per gene) plus a
`fitness(genes) -> float` function (lower is better). Everything
physics-specific — spectra, models, output files — stays behind that
interface; `Solvers` never reaches into it.

```python
from Solvers import get_solver
from Solvers.core import OptimizationProblem, ParameterSpace, GeneRange
import numpy as np

class ToyProblem(OptimizationProblem):
    def __init__(self, target):
        space = ParameterSpace([GeneRange.from_bounds(-10, 10, 0.01) for _ in target])
        super().__init__(space)
        self.target = np.asarray(target)

    def fitness(self, genes):
        return float(np.sum((genes - self.target) ** 2))

problem = ToyProblem(target=[1.0, 2.0, 3.0])
solver = get_solver("GA")(problem, options={"nPops": 100, "nGen": 100})
result = solver.run()
print(result.best_individual.genes, result.best_value)
```

`get_solver` looks up a class by name (case-insensitive) or by the numeric
`solOpt`/`solver_type` IDs EXAFS's `.ini` files have historically used
(`0` = GA, `1` = GA_Rechenberg, `2` = DE; `PSO` and `DE_MCMC` are
string-keyed only, no numeric ID — neither has a historical `solOpt`
value). `Solvers/__init__.py:SOLVER_REGISTRY` is the single source of
truth for both.

Optional hooks on `OptimizationProblem` — both no-ops by default:

- `sample_genes()` — override to bias initialization away from a uniform
  draw over the parameter space (e.g. DE-BNN's He-normal init).
- `on_generation_end(state, population)` / `on_run_end(state, population)` —
  problem-specific per-generation work (write output files, mid-run
  parameter sweeps — e.g. EXAFS's E0 optimization at `nGen // 2`).

## Layout

    core/
    ├── parameter_space.py   # GeneRange (discrete grid), ContinuousGeneRange
    │                        #   (continuous, for unbounded values like NN weights),
    │                        #   ParameterSpace (ordered collection of genes)
    ├── individual.py        # Individual: one genome (flat gene vector) + mutate()
    ├── population.py        # Population: generate/evaluate/sort a generation
    ├── problem.py            # OptimizationProblem ABC — the plug-in contract
    ├── run_state.py          # RunState: generation counter, timing, best-fit tracking
    ├── result.py              # SolverResult: per-generation history + best individual
    └── base_solver.py        # BaseSolver ABC: owns the generation loop, subclasses implement step()
    ga/                        # Genetic Algorithm
    ├── selectors.py          # RouletteWheelSelector (0), TournamentSelector (1, stub)
    ├── crossovers.py          # Uniform (0), SinglePoint (1), DualPoint (2, stub),
    │                        #   Arithmetic (3), Or (4), Average (5)
    ├── mutators.py            # PerIndividual (0), PerGene (1), PerTrait (2, stub),
    │                        #   Metropolis (3), Bounded (4)
    └── ga_solver.py          # GASolver, GARechenbergSolver (adds the 1/5 success rule)
    de/
    └── de_solver.py          # DESolver + differential_evolution_step(population, F, CR) —
                             #   a standalone function physics modules that run their own
                             #   generation loop (e.g. EXAFS) can call directly
    pso/
    └── pso_solver.py          # PSOSolver — global-best Particle Swarm Optimization
    demcmc/                    # DE-BNN: DE reinterpreted as an MCMC sampler for BNN training
    ├── mlp.py                # MLPStructure — flat-genome feedforward MLP
    ├── bnn_problem.py        # BNNRegressionProblem — reference OptimizationProblem
    ├── mutation.py            # rand/1, rand/2, best/1, best/2 + MCMC noise term
    ├── hyper_mutation.py      # StagnationTracker + sample_hyperparameters (adapts F/CR/operator)
    ├── refinement.py          # svd_refine, local_search_refine, cluster_refine (stub)
    ├── posterior.py           # PosteriorResult: posterior sampling + predict/credible_interval
    └── demcmc_solver.py       # DEMCMCSolver ties the above into one BaseSolver
    tests/                     # unittest, numpy-only — runs without any physics module installed

## Solvers at a glance

| name (`get_solver(...)`) | numeric ID | class |
|---|---|---|
| `GA` | 0 | `GASolver` |
| `GA_RECHENBERG` | 1 | `GARechenbergSolver` |
| `DE` | 2 | `DESolver` |
| `PSO` | — | `PSOSolver` |
| `DE_MCMC` | — | `DEMCMCSolver` |

- **GA**: selection → crossover → mutation → evaluate, each stage a
  swappable operator selected by numeric option ID (see `ga/` above), so
  existing physics-module `.ini` files (`selOpt`, `croOpt`, `mutOpt`) behave
  identically after routing through `Solvers`.
- **GA_Rechenberg**: `GASolver` plus `rechenberg_update` — adapts the
  mutation chance each generation based on how often recent generations
  improved the best fit (1/5 success rule).
- **DE**: DE/rand/1/bin (Storn & Price 1997) — for each target, build a
  trial from three distinct donors (`mutant = a + F*(b - c)`), binomial
  crossover with the target, greedy one-to-one replacement. Also exposed as
  the free function `differential_evolution_step(population, F, CR)` for
  physics modules that own their generation loop and just want the DE math
  per-generation (see `PhysicsModules/EXAFS/exafs_neo/neoSolver.py`).
- **PSO**: global-best Particle Swarm Optimization (Kennedy & Eberhart
  1995, Clerc & Kennedy 2002 constriction coefficients) — each individual
  carries a velocity pulled toward its own best-known position and the
  swarm's best-known position (`v = w*v + c1*r1*(pBest - x) + c2*r2*(gBest
  - x)`), then moves and is re-clipped into the parameter space. No
  historical `solOpt` precedent (string-keyed only, like `DE_MCMC`).
- **DE_MCMC (DE-BNN)**: differential evolution reinterpreted as an MCMC
  sampler for Bayesian neural network posterior sampling (Forbes & Long,
  *Neurocomputing* 678 (2026) 133103). Not tied to any `PhysicsModules`
  entry — a general-purpose regression/BNN capability that lives in
  `Solvers` per the paper's own future-work suggestion. See the root
  `README.md`'s "DE-BNN: Bayesian neural networks via DE-MCMC" section for
  the full write-up and a worked example, and `design.md` for what's
  implemented vs. the one documented stub (`cluster_refine`).

## Adding a new solver

Subclass `BaseSolver`, implement `step()` (advance the population by one
generation in place), and register it in `SOLVER_REGISTRY` in
`Solvers/__init__.py`. `BaseSolver.run()` owns the generation loop: it
initializes the population, then per generation calls `step()`, collects
into a `SolverResult`, and fires the problem's `on_generation_end` hook —
finishing with `on_run_end`. Physics modules pick up the new solver by name
with no code changes on their side.

## Tests

From the repository root (no physics dependencies needed):

    python -m unittest discover -s Solvers/tests -t . -v
