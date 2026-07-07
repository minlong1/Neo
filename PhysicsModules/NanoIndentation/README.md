# NanoIndentation Neo

Fits the unloading segment of nanoindentation load–displacement curves with
the NEO solvers, using the Oliver-Pharr power law:

    y = A (h - h_f)^m

Ported from the standalone `nano_indent` package (Nano_Neo); the genetic
algorithm that used to be embedded here now comes from the shared,
physics-agnostic `Solvers` package.

## Layout

    nanoindentation_neo/
    ├── nano_neo_data.py   # data loading (Hysitron / iMicro / plain CSV) + unloading-segment slicing
    ├── pathObj.py         # fit models (OliverPharr) and default parameter grids
    ├── problem.py         # NanoIndentationProblem — the Solvers plug-in point
    ├── parser.py          # .ini parsing/validation (same file format as nano_indent)
    ├── nano_indent.py     # NanoNeo runner: read -> setup -> run, logging + output files
    ├── input_arg.py       # CLI entry point (`nano_neo`)
    └── gui/               # legacy tkinter GUI, copied as-is (run scripts directly; not wired up)

## Usage

From the repository root:

    nano_neo -i my_run.ini

The `.ini` format is unchanged from nano_indent — sections `[Inputs]`
(`data_file`, `output_file`, `data_cutoff`), `[Populations]` (`population`,
`num_gen`, `best_sample`, `lucky_few`), `[Mutations]` (`chance_of_mutation`,
`mutated_options`), `[Paths]` (`npaths`, `fits`, optional
`a_range`/`hf_range`/`m_range` as `low, high, step`), `[Outputs]`
(`print_graph`, `num_output_paths`).

`mutated_options` maps onto the shared solvers as:

| value | behavior |
|-------|----------|
| 0     | Rechenberg-adaptive GA, mutation regenerates whole individuals (original default) |
| 1     | GA, per-gene mutation |
| 2     | GA, Metropolis mutation |
| 3     | GA, bounded mutation |

Outputs (per generation) go to `output_file` (fit history), `*_data.csv`
(best-fit parameters), and `*.log` (run log).

Programmatic use:

```python
from Solvers import get_solver
from PhysicsModules.NanoIndentation.nanoindentation_neo.problem import NanoIndentationProblem

problem = NanoIndentationProblem("curve.csv", data_cutoff=(0.1, 0.9), npaths=1)
result = get_solver("GA")(problem, options={"nPops": 500, "nGen": 100}).run()
print(result.best_individual.genes)  # [A, h_f, m]
```

## Tests

From the repository root (numpy-only; uses synthetic curves, no instrument
data needed):

    python -m unittest discover -s PhysicsModules/NanoIndentation/tests -t . -v

## Differences from nano_indent

- The GA (selection, crossover, mutation, Rechenberg rule) lives in
  `Solvers/`; this module only defines the physics (data + model + fitness).
- Genomes scoring NaN (h < h_f regions) are scored `inf` instead of being
  discarded and regenerated — same effect through selection, less machinery.
- Input parsing happens in functions at run time, not at module import time.
- The Rechenberg mutation-chance step uses the shared Solvers schedule
  (±0.25 percentage points per generation vs. the original ±0.5).
