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
    ├── model.py        # ModelSpec (expression + free/fixed/linked params by name), DEFAULT_MODEL_SPEC, build_xspec_model()
    ├── problem.py       # AstroNeoProblem — the Solvers plug-in point (loads the spectrum, scores genomes)
    ├── parser.py        # .ini parsing/validation
    ├── astro_neo.py      # AstroNeo runner: read -> setup -> run, logging + output files
    ├── input_arg.py      # CLI entry point (`astro_neo`)
    └── helper.py         # logger + banner (re-exports PhysicsModules.common)

## The default model, and fitting a different one

`AstroNeoProblem` fits one `model.ModelSpec` at a time — an XSPEC model
expression plus its free parameters (the fit genome), fixed/frozen
parameter overrides, and any parameter links, all addressed by dotted
`component.parameter` name (e.g. `"powerlaw.PhoIndex"`) rather than XSPEC's
raw positional parameter index. `model.DEFAULT_MODEL_SPEC` — used unless
you pass a different one — is the original ported model: an absorbed
powerlaw plus a two-temperature APEC plasma and an ACX2 charge-exchange
component,

    TBabs(TBabs*powerlaw + lsmooth(vapec + vapec + zashift*vacx2))

same expression, fixed values, links, and free-parameter set as the
original script, just re-expressed by name.

To fit a different model, construct a `ModelSpec` and pass it to
`AstroNeoProblem(model_spec=...)`:

```python
from PhysicsModules.AstroNeo.astro_neo.model import ModelSpec
from PhysicsModules.AstroNeo.astro_neo.problem import AstroNeoProblem

custom_spec = ModelSpec(
    expr="TBabs*powerlaw",
    free_params=("powerlaw.PhoIndex", "powerlaw.norm"),
    param_ranges={
        "powerlaw.PhoIndex": (0.5, 3.0),
        "powerlaw.norm": (1e-5, 1e-2),
    },
    fixed_params={"TBabs.nH": 0.0279},
    # unfrozen=(...,)  # dotted names to explicitly .frozen = False
    # links=(("target.param", "source.param"), ...)
)
problem = AstroNeoProblem(data_dir="...", data_file="spectrum.fits", model_spec=custom_spec)
```

Names, not indices, because indices are only meaningful for one specific
expression (renumber if it changes at all) and can otherwise only be
discovered by building the model and reading `model.show()` — fragile and
opaque for defining a new spec. Every `xspec.Parameter` exposes its own
`.index`, so `AstroNeoProblem.__init__` resolves `free_params` to indices
once, up front; `fitness()` still uses `Model.setPars({index: value})`,
XSPEC's fast bulk setter, every evaluation — resolving by name on *every*
call instead measured ~14x slower in practice (bulk index-based setPars:
~0.4s/eval; per-parameter name-based attribute sets: ~5.5s/eval), which
would undo the point of building the `Model` once per problem instance
(see "Differences from the original script" below). See `model.py`'s
module docstring for the full mechanics.

This is scoped to the Python API, not the `.ini` — a nested spec like this
doesn't map cleanly onto flat `.ini` sections, and the `astro_neo` CLI's
one job is running `DEFAULT_MODEL_SPEC` against a spectrum. A custom
model is a notebook/script use case for now.

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
optional `bg_file`/`rsp_file` — override the background/response `xspec.Spectrum()`
loads; leave unset to use whatever `data_file`'s own FITS header
(`BACKFILE`/`RESPFILE` keywords) already points at, which is PyXspec's
default and is correct for a normally-grouped PHA file — `acx2_path`,
`xmin`/`xmax` ignore-range bounds in Angstrom), `[Populations]` (`population`, `num_gen`, optional
`best_sample`/`lucky_few`), `[Mutations]` (`solver_type` — 0 GA, 1 GA with
Rechenberg, 2 DE, default DE — plus GA's `mutated_options`/
`crossover_options`/`chance_of_mutation` or DE's `mutf`/`mutcr`),
`[Outputs]` (`print_graph`).

Outputs (per generation) go to `output_file` (fit history), `*_data.csv`
(best-fit parameters by name), and `*.log` (run log). With `print_graph =
True`, `run()` plots the global-best individual against the data at the
end of the run (`AstroNeo.plot()`) — data with error bars vs. the folded
model, over PyXspec's own `Plot("data")` arrays and axis labels (matplotlib,
requires a display; not run under CI/headless verification).

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

- `DEFAULT_MODEL_SPEC`'s fixed/frozen values were tuned against one
  specific instrument/source; fitting a different target with the same
  model expression still means constructing a new `ModelSpec` with new
  fixed values (see "The default model, and fitting a different one").
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
- The model is no longer a single hardcoded expression — `model.py`'s
  `ModelSpec` generalizes it to any XSPEC expression, with parameters
  addressed by dotted name instead of the original's raw positional
  indices (see "The default model, and fitting a different one").
  `DEFAULT_MODEL_SPEC` reproduces the original model exactly — same
  expression, fixed values, links, and free-parameter set, verified to
  give the identical `Fit.statistic` — just re-expressed by name.
