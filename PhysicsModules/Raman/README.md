# Raman Neo

Fits Raman spectra with the NEO solvers as a sum of analytic spectral bands
(Lorentzian / Gaussian / pseudo-Voigt), optimising each band's amplitude,
width, and centre offset plus one global wavenumber shift.

Ported from the standalone [`raman_neo`](https://github.com/laumiulun/raman_neo)
package, which was itself a Raman re-implementation of EXAFS Neo. The
genetic algorithm that used to be embedded there (a fork of the EXAFS one)
now comes from the shared, physics-agnostic `Solvers` package, the same way
NanoIndentation's was replaced.

## Relationship to EXAFS

The physics maps one-to-one onto the EXAFS module, so the genome has the
same shape: one shared gene plus three genes per band/path.

| EXAFS                         | Raman                              |
| ----------------------------- | ---------------------------------- |
| χ(k) vs wavenumber k          | Intensity vs Raman shift (cm⁻¹)    |
| FEFF scattering path          | Spectral band (line shape)         |
| S0² (amplitude)               | `amp`: band amplitude              |
| E0 (shared energy shift)      | `shift`: shared wavenumber offset  |
| σ² (Debye–Waller)             | `width`: band width (cm⁻¹)         |
| ΔR (radial shift)             | `center_shift`: per-band offset    |
| Larch/FEFF `path2chi`         | analytic line-shape synthesis      |
| k-weighting                   | optional wavenumber weighting      |
| `autobk` background           | polynomial baseline subtraction    |
| path optimize (area)          | band optimize (integrated area)    |

Each band's absolute centre is `nominal_center + center_shift + shift`.
Genome layout (flat, see `problem.py`):

    genes = [shift, amp_0, width_0, cshift_0, amp_1, width_1, cshift_1, ...]

Like EXAFS's E0, the shared shift is swept over its whole grid at the run's
midpoint and at the end (`RamanProblem.on_generation_end`/`on_run_end`), and
after the midpoint sweep, newly sampled genomes start from the best shift
found. Set `shiftSweep: False` in the parameter dict (or
`RamanProblem(shift_sweep=False)`) to turn that off.

The model is synthesised directly with NumPy, so the core fit needs
**numpy only**. Matplotlib is needed only for `print_graph = True` /
`RamanAnalysis.plot()`.

## Layout

    raman_neo/
    ├── bandshape.py       # line shapes, Band, band-definition file reader, polynomial baseline
    ├── bandrange.py       # per-band (amp, width, center_shift) grids + shared-shift grid, custom range file
    ├── raman_data.py      # spectrum loading, resampling onto the model grid, baseline, normalisation
    ├── problem.py         # RamanProblem: the Solvers plug-in point (fitness + shift-sweep hooks)
    ├── parser.py          # .ini parsing/validation (same file format as raman_neo)
    ├── raman.py           # RamanNeo runner: read -> setup -> run, logging + output files
    ├── input_arg.py       # CLI entry point (`raman_neo`)
    ├── helper.py          # logger + banner (re-exports PhysicsModules.common)
    ├── analysis/
    │   ├── raman_score.py     # standalone synthesis/scoring (larch_score analogue)
    │   └── raman_analysis.py  # RamanAnalysis: peak table, R-factor, report, plot
    └── data/              # synthetic sample spectrum, band file, example .ini

## Usage

Install the module's extra (matplotlib, only for plotting) from the
repository root:

    pip install -e ".[raman]"

Run the bundled example (synthetic spectrum with bands at 480/720/1350/1580
cm⁻¹; the band file has two decoy rows that `band_optimize` should flag):

    raman_neo -i PhysicsModules/Raman/raman_neo/data/sample_raman.ini

### Input `.ini`

Unchanged from raman_neo. Sections:

- `[Inputs]`: `data_file`, `output_file`, `band_file`; optional
  `num_compounds` (default 1; with more than one, `band_file` is a
  comma-separated list, one file per compound) and `bandrange_file`.
- `[Populations]`: `population`, `num_gen`, `best_sample`, `lucky_few`.
  `best_sample`/`lucky_few` are **counts** out of `population` (raman_neo's
  convention), not percents.
- `[Solver]`: optional `solver_options`: `0` GA (default), `1` GA +
  Rechenberg, `2` DE, or any `Solvers` registry name (`PSO`, `CMA_ES`).
  Optional `mutf`/`mutcr` set DE's F/CR (default 0.5/0.9).
- `[Mutations]`: `chance_of_mutation` (percent); optional
  `mutated_options` (default 1), `selection_options`, `crossover_options`.
- `[Bands]`: `band_list`: 1-based rows of the band file to fit (`1,2,3`, or
  `[1,2][3,4]` per compound); optional `band_optimize`,
  `band_optimize_percent` (default 0.01).
- `[Spectrum]`: `xmin`, `xmax` (fit window, cm⁻¹), `weight` (residuals ×
  x^weight; 0 = unweighted), `deltax` (model grid step); optional
  `baseline_order` (0 = none), `normalize` (default True: unit peak in the
  window), `shift_range` (`low, high, step`, default `-10, 10.05, 0.1`).
- `[Outputs]`: `print_graph`; optional `steady_state_exit` (parsed, unused,
  as in NanoIndentation).

`original_chance_of_mutation`, `chance_of_mutation_e0`, `individual_band`,
`optimize_only`, and `num_output_bands` are still accepted so existing
raman_neo `.ini` files run unchanged, but they have no effect (see
"Differences from raman_neo").

The option IDs map straight onto `Solvers`:

| key                 | value | behaviour |
|---------------------|-------|-----------|
| `solver_options`    | 0 / 1 / 2 | GA / GA + Rechenberg / DE (`get_solver` numeric IDs) |
| `mutated_options`   | 0 | per individual (regenerate whole genomes) |
|                     | 1 | per gene (raman_neo's "per band") |
|                     | 2 | per trait (stub, no-op) |
|                     | 3 | Metropolis |
|                     | 4 | bounded |
| `selection_options` | 0 | roulette wheel (1 = tournament, stub) |
| `crossover_options` | 0–5 | uniform, single point, dual point (stub), arithmetic, or, average |

### Data file

Two columns, Raman shift (cm⁻¹) and intensity, whitespace- or
comma-separated. Lines starting with `#` are comments; rows needn't be
sorted.

### Band-definition file

One candidate band per line (the analogue of a FEFF path folder). Only the
centre is required:

```
# center, shape, eta, base_width
480.0,  lorentzian, 1.0, 9.0
1350.0, voigt,      0.7, 28.0
```

`shape` is `lorentzian`, `gaussian`, or `voigt`. `eta` is the Lorentzian
fraction of a pseudo-Voigt. `base_width` is a reference width in cm⁻¹
(informational). `width` is the Lorentzian HWHM or the Gaussian σ.

### Custom ranges

`bandrange_file`: one row per fitted band, 10 comma-separated columns:
`band, amp_low, amp_high, amp_dt, width_low, width_high, width_dt,
cshift_low, cshift_high, cshift_dt`. Defaults are amp 0–1.98 (step 0.02,
plus 0 so a band can switch off), width 2–39 cm⁻¹ (step 1), and
center_shift −15 to 15 cm⁻¹ (step 1).

### Outputs

Next to `output_file`:

- `output_file`: per-generation fit history (`Gen,TPS,FITTNESS,CURRFIT,CURRIND,BESTFIT,BESTIND`).
- `*_data.csv`: per generation, the best fit's `[amp, shift, width, center_shift]` per band (EXAFS layout).
- `*_peaks.csv`: final peak table (fitted centre, amplitude, width, area, area fraction per band).
- `*.log`: the run log, ending with the fit report.

### Programmatic use

```python
from Solvers import get_solver
from PhysicsModules.Raman.raman_neo.problem import RamanProblem
from PhysicsModules.Raman.raman_neo.analysis.raman_analysis import RamanAnalysis

problem = RamanProblem(
    "spectrum.dat", "bands.txt", band_list=[1, 2, 3, 4],
    xmin=150, xmax=1750, baseline_order=1,
)
result = get_solver("GA_Rechenberg")(problem, options={"nPops": 120, "nGen": 40}).run()

analysis = RamanAnalysis(problem, result.best_individual.genes)
print(analysis.report())
analysis.write_peak_table("peaks.csv")
```

The final shift sweep runs in `on_run_end`, after `BaseSolver.run()` has
collected the last generation, so it updates `result.best_individual.genes`
(the same object) and `solver.state.globBestVal`, but not
`result.best_value`. `RamanAnalysis` recomputes the fitness from the genes
unless you pass one. The `RamanNeo` runner re-syncs `result.best_value`
itself.

Or drive the `.ini` workflow from Python:

```python
from PhysicsModules.Raman.raman_neo.raman import RamanNeo

raman = RamanNeo()
raman.raman_read(filepath="my_run.ini")   # or input_parameters={...}
raman.raman_setup()
result = raman.run()
print(raman.analysis.report())
```

## Tests

From the repository root (numpy only; synthetic spectra, no instrument data
needed):

    python -m unittest discover -s PhysicsModules/Raman/tests -t . -v

## Differences from raman_neo

- The GA (selection, crossover, mutation, Rechenberg rule) lives in
  `Solvers/`. The module only defines the physics (data + model + fitness).
  raman_neo's `neoPars`/`neoSelector`/`neoCrossOver`/`neoMutator`/
  `neoSolver`/`neoResult`/`individual`/`bandObj` are gone. Nothing depended on
  their internals (raman_neo was at v0.1.0), so there's no back-compat shim
  layer like EXAFS's.
- `solver_options = 2` now runs real differential evolution
  (`Solvers.de`). In raman_neo it was a no-op stub, as it was in upstream
  EXAFS Neo. PSO and CMA-ES are also available by name.
- The shared shift is an ordinary gene: the selected mutator and
  crossover act on it like any other. raman_neo's separate
  second-half-only shift mutation (`chance_of_mutation_e0`) is dropped. The
  mid-run and final shift sweeps are kept.
- The shift sweep covers the shift gene's own grid (default ±10 cm⁻¹). In
  raman_neo it scanned ±30 cm⁻¹ while the gene was only ever sampled from
  ±10. Widen `shift_range` if your calibration offset can be larger.
- The polynomial baseline is fitted only over grid points the spectrum
  actually covers. raman_neo fitted it over the whole 0..xmax model grid,
  where everything below the first data point is zero padding, which
  dragged the lower-envelope fit down.
- Input parsing happens in functions at run time, not through an
  `InputParamsParser` object.
- New outputs: `*_peaks.csv`, and a fit report at the end of the log.
- Genomes scoring NaN are scored `inf`, as in NanoIndentation.
