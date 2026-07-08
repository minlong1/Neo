# NEO

NEO is a modular scientific computing framework. Optimization algorithms live
in a shared, physics-agnostic `Solvers` package; each physics application
lives in its own module under `PhysicsModules` and plugs into the solvers by
implementing a small problem interface.

## Repository Structure

    NEO/
    ├── Solvers/                  # physics-agnostic solvers (numpy only)
    │   ├── core/                 # ParameterSpace, Individual, Population,
    │   │                         #   OptimizationProblem, BaseSolver, ...
    │   ├── ga/                   # Genetic Algorithm (+ Rechenberg variant)
    │   ├── de/                   # Differential Evolution (DE/rand/1/bin)
    │   ├── demcmc/                # DE-BNN: DE as an MCMC sampler for Bayesian NNs
    │   └── tests/
    ├── PhysicsModules/
    │   ├── EXAFS/                # EXAFS Neo — fully supported
    │   ├── NanoIndentation/      # Nano Neo — fully supported
    │   └── XPS/                  # XPS Neo — fully supported (self-contained GA, see its README)
    ├── pyproject.toml
    └── README.md

## Physics modules

**EXAFS** fits Extended X-ray Absorption Fine Structure data with a genetic
algorithm. Setup instructions, dependencies, example commands, and test
information are in
[PhysicsModules/EXAFS/README.md](PhysicsModules/EXAFS/README.md).

**NanoIndentation** (Nano Neo) fits the unloading segment of nanoindentation
load–displacement curves with the Oliver-Pharr model. See
[PhysicsModules/NanoIndentation/README.md](PhysicsModules/NanoIndentation/README.md).

**XPS** (XPS Neo) fits X-ray Photoelectron Spectroscopy spectra (Voigt,
Gaussian, Lorentzian, Double Lorentzian, Doniach-Sunjic peaks over Baseline/
Linear/Shirley/SVSC/Tougaard backgrounds). Unlike the other two modules, it
keeps its own GA/DE loop rather than routing through `Solvers` — its genome
is a heterogeneous, type-tagged list, not a fixed-width float vector; see
[PhysicsModules/XPS/README.md](PhysicsModules/XPS/README.md) for why and for
its (pytest-based) test suite.

## Installing a single module

Each physics module's third-party dependencies are an optional extra, so
installing one doesn't drag in the others' (NanoIndentation shouldn't need
`xraylarch`, and neither should need `numba`). The bare install (no extras)
gets you `Solvers` alone — no physics, `numpy` only.

From the repository root:

```bash
# Solvers only (e.g. developing/using the optimization framework directly)
pip install -e .

# EXAFS: adds xraylarch, scipy, attrs, matplotlib, psutil
pip install -e ".[exafs]"
exafs_neo -i <your_exafs_input.ini>
exafs_neo_gui

# NanoIndentation: adds matplotlib (only needed for printGraph = True;
# the core fit runs on numpy alone)
pip install -e ".[nanoindentation]"
nano_neo -i <your_nanoindent_input.ini>

# XPS: adds scipy, numba, matplotlib, psutil
pip install -e ".[xps]"
xps_neo -i <your_xps_input.ini>
xps_neo_gui

# Everything (all three modules' extras at once)
pip install -e ".[all]"
```

Extras are additive — installing a second module's extra doesn't uninstall
the first's. Each module's own README has more detail (input file format,
GUI notes, test commands): `PhysicsModules/EXAFS/README.md`,
`PhysicsModules/NanoIndentation/README.md`, `PhysicsModules/XPS/README.md`.

## The Solvers contract

A physics module exposes its fitting task by subclassing
`Solvers.core.OptimizationProblem` — a `ParameterSpace` describing the fit
parameters plus a `fitness(genes)` function — and can then run any registered
solver:

```python
from Solvers import get_solver

solver = get_solver("GA")(problem, options={"nPops": 100, "nGen": 100})
result = solver.run()
```

Available solvers: `GA`, `GA_Rechenberg`, `DE`, `DE_MCMC`. New solvers are
added under `Solvers/` and registered in `Solvers/__init__.py` — physics
modules pick them up without code changes. `Solvers.de` also exposes
`differential_evolution_step(population, F, CR)` as a standalone function,
for physics modules (like EXAFS) that run their own generation loop and
call it per-generation rather than handing control to a `BaseSolver`.

### DE-BNN: Bayesian neural networks via DE-MCMC

`Solvers.demcmc` implements DE-BNN (Forbes & Long, "DE-BNN: An
evolutionary approach to Bayesian neural network posterior sampling",
*Neurocomputing* 678 (2026) 133103): differential evolution reinterpreted
as an MCMC sampler. DE's own selection rule (accept a trial if it scores
better) already matches an MCMC acceptance step; with small Gaussian noise
added to the mutant for detailed balance, the sequence of accepted
candidates at each population index becomes a Markov chain, giving up to
`nPops` chains "for free". Applied to training a neural network's weights
and biases, the population of post-burn-in samples *is* the posterior
distribution — used here for Bayesian regression with a mean prediction,
credible interval, and a mode point estimate, rather than a single trained
network.

```python
from Solvers.demcmc import BNNRegressionProblem, DEMCMCSolver

problem = BNNRegressionProblem(layer_sizes=[7, 35, 10, 5, 1], X=X_train, y=y_train)
solver = DEMCMCSolver(problem, options={
    "nPops": 200, "nGen": 12000, "burn_in": 7000, "num_chains": 4,
    "F": 0.5, "CR": 0.7, "mutation_operator": "rand/1",
    "hyper_mutation": True, "svd": True, "local_search": True,
})
solver.run()

mean_prediction, _ = solver.posterior.predict(problem, X_test)
lo, hi = solver.posterior.credible_interval(problem, X_test, low=5, high=95)
```

Implemented from the paper: the base DE-MCMC loop with configurable
mutation operator (`rand/1`, `rand/2`, `best/1`, `best/2`) and MCMC noise
term (Section 2, Eq. 48); hyper-mutation, which resamples F/CR/operator
when a running-average residual signals stagnation (Section 3.1-3.4); SVD
and local-search refinement, run periodically and accepted only if they
beat the population's current best (Section 3.5, 3.7); single- and
multi-chain posterior collection and prediction (Section 5). **Not**
implemented: clustering-based refinement (Section 3.6) — it would need a
clustering dependency (e.g. scikit-learn) for one of three interchangeable
refinement techniques; `Solvers.demcmc.cluster_refine` is a documented
stub raising `NotImplementedError`.

## Running tests

From the repository root:

    # solver framework (no physics dependencies needed)
    python -m unittest discover -s Solvers/tests -t . -v

    # EXAFS (requires xraylarch; see the EXAFS README)
    cd PhysicsModules/EXAFS && python -m unittest discover -s tests -t ../.. -v

    # NanoIndentation (numpy only)
    python -m unittest discover -s PhysicsModules/NanoIndentation/tests -t . -v

    # XPS (pytest, not unittest — see PhysicsModules/XPS/README.md)
    cd PhysicsModules/XPS && pytest -m "not golden" -q

Typical runs from the repository root:

    exafs_neo -i <your_exafs_input.ini>
    nano_neo -i <your_nanoindent_input.ini>
    xps_neo -i <your_xps_input.ini>
