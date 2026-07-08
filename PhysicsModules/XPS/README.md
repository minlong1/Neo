# XPS Neo

Genetic-algorithm curve fitting for X-ray Photoelectron Spectroscopy (XPS)
spectra. Fits peak models (Voigt, Gaussian, Lorentzian, Double Lorentzian,
Doniach-Sunjic; singlets, spin-orbit doublets, Coster-Kronig) over physical
backgrounds (Baseline, Linear, Shirley, SVSC, Tougaard, ...) using a
population-based optimizer (classic GA or differential evolution).

Ported from the standalone `XPS_Neo` package, which had already gone through
five phases of production-readiness work (see `design.md` in that repo):
determinism/seeding, a golden-master test harness, killing import-time side
effects, unifying a GUI/CLI fork, and a performance pass — all with bit-exact
regression discipline. This port carries that work over essentially
unchanged; see "Why this module doesn't use `Solvers`" below for the one
structural decision that differs from EXAFS/NanoIndentation.

## Usage

```bash
xps_neo -i <config.ini> [--seed N] [--workers N]
```

- `--seed N` — reproducible runs (seeds all RNGs).
- `--workers N` — parallel population evaluation; results are identical to
  serial.
- `--version`, `-v` (verbose), `-s` (echo parsed input).

The INI format is documented by example in `tests/golden/cases/*.ini`
(working, tested configurations — the old sample INIs under
`xps_neo/gui/*.ini` predate required keys and no longer run). Each run
writes `<output>.txt` (per-generation log), `<output>_data.csv` (best-fit
parameters per generation), and `<output>.log`.

GUI: `xps_neo_gui` (requires tkinter). It builds INI files from form fields
and launches the installed `xps_neo` CLI as a subprocess — if you edit
`xps_neo/` during development, reinstall with `pip install -e ".[xps]"` from
the repo root or the GUI will keep running the old code.

## Layout

    xps_neo/
    ├── xps_data.py        # spectrum loader (Hysitron-style txt, offsets)
    ├── xps_fit.py         # peak/background shape math (numba-jitted)
    ├── xps_individual.py  # Individual: one candidate fit (peaks + backgrounds)
    ├── loss.py            # the one loss function (compute_loss), shared by
    │                       #   the GA and the GUI's post-analysis
    ├── parser2.py          # .ini -> section dict (named to avoid shadowing stdlib parser)
    ├── ini_parser.py       # section dict -> flat config values (load_config)
    ├── periodic_table.py / periodic_table_data.py  # literature BE/width lookup for the GUI
    ├── xps.py              # XPS_GA: GA/DE loop, fitness, peak add/remove, CLI main()
    ├── input_arg.py        # argparse + seeding
    └── gui/                # xps_neo_gui.py, xps_analysis2.py, xps_plot.py (unified GUI)

## Why this module doesn't use `Solvers`

EXAFS and NanoIndentation both fit a **flat, fixed-width vector of floats**
per generation, so their `Individual` is a thin wrapper around
`Solvers.core.Individual`/`ParameterSpace`, and their selection/crossover/
mutation come from `Solvers.ga`.

XPS's genome doesn't fit that shape. `Individual.get_params()` returns a
list of floats **interleaved with type-name strings** (`'Voigt'`,
`'Baseline'`, ...) that mark where one peak's or background's parameter
block ends — the number of floats per block depends on the peak type
(singlet vs. doublet vs. Coster-Kronig adds fields), and correlated/limited
parameters change which fields actually vary. `xps.py`'s own crossover
walks this heterogeneous list by finding the string markers
(`XPS_GA.crossover`'s `dividers`). Forcing this through `Solvers.core`'s
homogeneous-float-array model would mean either flattening peak-type
identity out of the genome into a separate fixed side-channel (a real
redesign of `xps_individual.py`/`xps_fit.py`, risking the bit-exact
behavior the golden-master suite pins) or building a second, parallel
heterogeneous-genome abstraction inside `Solvers` for one caller.

Given `Solvers.ga` doesn't (yet) buy XPS anything it doesn't already have —
its own GA has Rechenberg, random-perturbation, and Metropolis mutation
options plus a working, jDE-adapted differential evolution solver, already
characterized by the golden matrix — the module keeps its own GA/DE loop
(`XPS_GA` in `xps.py`) rather than routing through `Solvers.core`. If a
future change needs a solver feature that lives only in `Solvers` (e.g. a
new crossover operator), that's the trigger to revisit this, not before.

## Tests

Unlike EXAFS/NanoIndentation (`unittest`, run via `python -m unittest
discover`), this module's tests are **pytest**-based — ported as-is along
with the code, rather than rewritten, so the bit-exact golden-master
discipline they encode carries over intact. From `PhysicsModules/XPS/`:

```bash
pip install -e "../..[xps]" -r requirements-dev.txt   # from repo root: pip install -e ".[xps]" -r PhysicsModules/XPS/requirements-dev.txt
pytest -m "not golden" -q     # unit + component tests (~2 s)
pytest -m golden -q           # golden-master matrix, runs the CLI per case (~1 min)
```

**Environment note:** the golden matrix is bit-exact only within the pinned
environment in `constraints.txt` (numpy 2.3.0/numba 0.63.1/scipy 1.17.0,
Python 3.14). In this repo's default environment (numpy 1.26.4/numba
0.59.0/scipy 1.15.3, Python 3.10), expect two kinds of environment-only
divergence, verified during the port to be display/precision artifacts, not
logic differences:

- 2 of 36 component tests fail at the ~1e-14 relative level (numba/numpy
  float-op reordering between versions).
- 5 of 9 golden cases fail: `assert text == expected` because numpy ≥2.0
  reprs scalars as `np.float64(1.23)` inside a printed list, where numpy
  1.26 prints plain `1.23` — plus, in the more numerically involved cases
  (Coster-Kronig, correlated peaks, Doniach-Sunjic/Shirley), the very last
  significant digit of some floats. Every case was manually re-verified
  with the `np.float64(...)` wrapper stripped and reproduced its expected
  trajectory generation-for-generation; see the port's commit message for
  the verification method. `lorentzian_baseline`, `voigt_baseline_svsc`,
  and `voigt_de_mode` are bit-exact even in this environment.

This is the environment-drift risk `design.md`'s own risk table anticipates
("numpy upgrade changes float behavior") — not a reason to touch the
comparator (design.md: "Never adjust the comparator to make a diff pass").
Run the golden matrix in the pinned `constraints.txt` environment for a
true pass/fail signal.

## Known limitations (carried over as-is)

Ported unchanged from the source project's own bug register — not fixed
here, since fixing them is an intentional, reviewed act with its own golden
regen, not a side effect of porting:

- Tougaard (2-/3-param) backgrounds crash (`known_broken/`).
- `peak_add_remove` (residual-driven peak addition) is unusable from the
  CLI (crashes on the first added peak).
- Only the background name `SVSC` is accepted; pre-rename
  `Shirley-Sherwood` INIs no longer work.

See `tests/golden/known_broken/README.md` for the quarantined cases.
