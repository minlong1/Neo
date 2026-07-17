# Astro Neo

Fits an X-ray CCD spectrum with the NEO solvers, using [PyXspec](https://heasarc.gsfc.nasa.gov/xanadu/xspec/python/html/)
to evaluate an absorbed powerlaw plus a two-temperature APEC plasma and an
[ACX2](https://github.com/AtomDB/ACX2) charge-exchange component:

    TBabs(TBabs*powerlaw + lsmooth(vapec + vapec + zashift*vacx2))

scored by the Cash statistic. Ported from a personal research script
(`Astro_Neo`); the genetic algorithm/DE that used to be embedded here now
comes from the shared, physics-agnostic `Solvers` package, the same as
EXAFS and NanoIndentation.

## Layout

    astro_neo/
    ├── model.py        # the one fitted model: expression, free-parameter bounds, fixed/linked params, build_xspec_model()
    ├── problem.py       # AstroNeoProblem — the Solvers plug-in point (loads the spectrum, scores genomes)
    ├── parser.py        # .ini parsing/validation
    ├── astro_neo.py      # AstroNeo runner: read -> setup -> run, logging + output files
    ├── input_arg.py      # CLI entry point (`astro_neo`)
    └── helper.py         # logger + banner (re-exports PhysicsModules.common)

## Why the model is hardcoded

Unlike EXAFS (any set of FEFF paths) or NanoIndentation (any number of
Oliver-Pharr terms), this module fits exactly one XSPEC model expression —
same as the original script. XSPEC's `Model.setPars(dict)` only accepts
1-based positional parameter indices, not names, so `model.py`'s
`FREE_PARAM_INDICES` and `FIXED_PARAMS` are tied to `MODEL_EXPR` exactly as
written; changing the expression means regenerating both (e.g. via
`model.show()` in a PyXspec session). Generalizing to an arbitrary
user-supplied XSPEC model string is a real feature, not attempted here —
if that's needed, it's the trigger to add it deliberately.

## Setup

PyXspec ships with a HEASOFT source build, not PyPI or conda-forge — there
is no `pip install -e ".[astro]"` for this module. You need:

1. A working HEASOFT/XSPEC install (built from source; see the
   [HEASoft docs](https://heasarc.gsfc.nasa.gov/docs/software/heasoft/)).
2. `HEADAS`/`ATOMDB` exported and `headas-init.sh` sourced *before* starting
   Python, so `import xspec` resolves — e.g.:

   ```bash
   export ATOMDB=$HOME/atomdb
   export HEADAS=/path/to/heasoft-<ver>/<platform>
   source "$HEADAS"/headas-init.sh
   ```

3. (Optional) [ACX2](https://github.com/AtomDB/ACX2)'s `acx2_xspec.py` on
   `sys.path` if the model needs the charge-exchange component — pass its
   directory as the `.ini` `acx2_path` key (or set it before constructing
   `AstroNeoProblem` directly); this is a separate, external plasma-physics
   model with its own AtomDB data dependency, not vendored into this repo
   (same treatment as EXAFS's optional `contrib/sabcor` submodule).

`problem.py`/`model.py` import `xspec`/`acx2_xspec` lazily (inside
`__init__`/`fitness`, not at module scope), so the rest of the module stays
importable — and its non-PyXspec tests runnable — without any of this set up.

## Usage

From the repository root, with the environment above active:

    astro_neo -i my_run.ini

`.ini` sections: `[Inputs]` (`data_dir`, `data_file`, `output_file`,
optional `bg_file`/`rsp_file` — accepted but not yet wired into the
spectrum load, matching upstream — `acx2_path`, `xmin`/`xmax` ignore-range
bounds in Angstrom), `[Populations]` (`population`, `num_gen`, optional
`best_sample`/`lucky_few`), `[Mutations]` (`solver_type` — 0 GA, 1 GA with
Rechenberg, 2 DE, default DE — plus GA's `mutated_options`/
`crossover_options`/`chance_of_mutation` or DE's `mutf`/`mutcr`),
`[Outputs]` (`print_graph`).

Outputs (per generation) go to `output_file` (fit history), `*_data.csv`
(best-fit parameters by name), and `*.log` (run log).

Programmatic use:

```python
from Solvers import get_solver
from PhysicsModules.AstroNeo.astro_neo.model import PARAM_NAMES
from PhysicsModules.AstroNeo.astro_neo.problem import AstroNeoProblem

problem = AstroNeoProblem(data_dir="...", data_file="spectrum.fits")
result = get_solver("DE")(problem, options={"nPops": 20, "nGen": 20}).run()
print(dict(zip(PARAM_NAMES, result.best_individual.genes)))
```

## Tests

From the repository root:

    python -m unittest discover -s PhysicsModules/AstroNeo/tests -t . -v

The parameter-space/genome-mapping tests are pure numpy and always run.
The end-to-end fitness test needs PyXspec plus a real spectrum — neither is
available in CI — so it's skipped unless `ASTRO_NEO_TEST_DATA_DIR` and
`ASTRO_NEO_TEST_DATA_FILE` (and optionally `ASTRO_NEO_ACX2_PATH`) point at
one, e.g. for a manual local run:

    ASTRO_NEO_TEST_DATA_DIR=/path/to/spectra ASTRO_NEO_TEST_DATA_FILE=spectrum.fits \
        python -m unittest PhysicsModules.AstroNeo.tests.test_problem.TestAstroNeoProblemLive -v

## Known limitations (carried over from the source script)

- The fitted model expression is fixed (see "Why the model is hardcoded").
- `bg_file`/`rsp_file` are accepted in the `.ini` and threaded through to
  `AstroNeoProblem`, but not yet passed to `xspec.Spectrum()` — the
  original script had the same gap.
- `model.py`'s `FIXED_PARAMS`/`FREE_PARAM_INDICES` were tuned against one
  specific instrument/source; reusing this module for a different target
  needs new fixed values and possibly a different free-parameter set.
- No bundled test spectrum ships with this repo (unlike EXAFS's committed
  `cu_test_files/`) — real `.fits` data plus a HEASOFT build are both
  personal-machine setup, so full runs can only be verified manually, not
  in CI. See "Tests" above.

## Differences from the original script

- The GA/DE (selection, crossover, mutation) lives in `Solvers/`; this
  module only defines the physics (model + spectrum + fitness). The
  original's own duplicate `NeoSelector`/`NeoCrossover`/`NeoMutator`/
  `NeoSolver`/`NeoPopulations` framework is gone — most of its GA
  crossover/mutation operators referenced an `Individual` API (`get_e0`,
  `get_path`, `exafsPathPars`) that doesn't exist on this module's own
  `Individual`, i.e. were never actually working; DE (`solver_type=2`) was
  the one path exercised by the original's own demo, and is this module's
  default.
- Dropped: the large dead shape-function zoo in the original `pathObj.py`
  (Gaussian/Voigt/Doniach/Shirley/... — leftover, unused overlap with the
  `XPS` module's peak shapes), the unused `sklearn`-based DE clustering
  solver variant, and the personal absolute paths (`/Users/.../Astro_Neo/...`)
  previously hardcoded into the demo and default configs — all data/ACX2
  paths are now `.ini`-configurable.
- The spectrum is loaded once per `AstroNeoProblem` instance instead of
  being reloaded every generation (the original's serial evaluation path
  called its per-worker spectrum-load function on every generation; its
  distributed/multiprocessing path loaded once per worker — this keeps the
  once-per-run behavior, since re-loading the same static spectrum
  mid-run can't change the fit).
- The XSPEC `Model` object is also now built once per `AstroNeoProblem`
  instance and reused across every `fitness()` call (`model.setPars(...)`
  per evaluation, not a fresh `xspec.Model(...)` each time) — unlike the
  spectrum-load fix above, this *does* diverge from the original, which
  rebuilt the model from scratch per individual. That was ~25s/evaluation,
  making any run at a practical population size (`nPops >= 100`)
  impractically slow (hours-to-days); reusing the Model is PyXspec's
  standard usage pattern (`Fit.statistic` is a live-recomputed property)
  and safe here because every free parameter is set on every call, so
  there's no stale-value carryover between generations.
