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
    │   ├── de/                   # Differential Evolution (stub)
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

Available solvers: `GA`, `GA_Rechenberg`, `DE` (stub). New solvers are added
under `Solvers/` and registered in `Solvers/__init__.py` — physics modules
pick them up without code changes.

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
