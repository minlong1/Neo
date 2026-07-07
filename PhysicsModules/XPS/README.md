# XPS Neo (scaffold)

Module for fitting X-ray Photoelectron Spectroscopy (XPS) data with the NEO
solvers. **Not implemented yet** — this directory holds the skeleton a future
implementation drops into.

## Module contract

Every NEO physics module plugs into the shared, physics-agnostic solvers in
`Solvers/` by implementing `Solvers.core.OptimizationProblem`:

1. **Problem** (`xps_neo/problem.py`) — subclass `OptimizationProblem`:
   - build a `ParameterSpace` (one `GeneRange` per fitting parameter),
   - load experimental data in `__init__` (all I/O stays in the module),
   - implement `fitness(genes) -> float` (lower is better),
   - optionally override `sample_genes()` and the `on_generation_end` /
     `on_run_end` hooks for initialization bias and per-generation outputs.
2. **Input parsing** — add a parser for the module's input format and a small
   CLI entry point, then register it in the root `pyproject.toml` under
   `[project.scripts]` (e.g. `xps_neo = "PhysicsModules.XPS.xps_neo.input_arg:main"`).
3. **Tests** (`tests/`) — keep module tests here; run from the repository
   root with:

        python -m unittest discover -s PhysicsModules/XPS/tests -t . -v

Running a solver against the problem:

```python
from Solvers import get_solver
from PhysicsModules.XPS.xps_neo.problem import XPSProblem

problem = XPSProblem()
solver = get_solver("GA")(problem, options={"nPops": 100, "nGen": 100})
result = solver.run()
```

See `PhysicsModules/EXAFS/` for the complete, working reference
implementation of this contract.
