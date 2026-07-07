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
    │   ├── XPS/                  # scaffold, not implemented yet
    │   └── NanoIndentation/      # scaffold, not implemented yet
    ├── pyproject.toml
    └── README.md

## Physics modules

**EXAFS** is the working module: it fits Extended X-ray Absorption Fine
Structure data with a genetic algorithm. Setup instructions, dependencies,
example commands, and test information are in
[PhysicsModules/EXAFS/README.md](PhysicsModules/EXAFS/README.md).

**XPS** and **NanoIndentation** are scaffolds: their directories document the
module contract and contain placeholder problem classes and tests, ready for
a future implementation. They are not supported yet.

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

A typical EXAFS run from the repository root:

    exafs_neo -i <your_input.ini>
