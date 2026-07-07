"""Canonicalize XPS_Neo run outputs for golden-master comparison (design.md §3.1).

A run in directory `outdir` with output_file = fit_out.txt produces:
  fit_out.txt       per-generation CSV log (contains a wall-clock TPS field)
  fit_out_data.csv  best-fit parameter rows per generation
  fit_out.log       logger output (timestamps, durations — not compared)

Canonical form drops every time-derived field and keeps everything else
verbatim, so comparison is bit-exact on the deterministic content.
"""

import pathlib

OUTPUT_STEM = "fit_out"


def canon_generations(text):
    """Drop the TPS (time-per-step) field from each generation line.

    Line format: gen,TPS,currfit,currind...),bestfit,bestind...
    Only the first two comma-separated fields are scalar; the parameter
    blobs after them contain commas, so split exactly twice.
    """
    out = []
    for line in text.splitlines():
        if line.startswith("Gen,"):  # header
            out.append(line)
            continue
        gen, _, rest = line.partition(",")
        _tps, _, payload = rest.partition(",")
        out.append(gen + "," + payload)
    return "\n".join(out) + "\n"


def canonicalize_run(outdir):
    """Return {relative filename: canonical text} for a completed run."""
    outdir = pathlib.Path(outdir)
    gen_file = outdir / f"{OUTPUT_STEM}.txt"
    data_file = outdir / f"{OUTPUT_STEM}_data.csv"
    return {
        "generations.txt": canon_generations(gen_file.read_text()),
        "params.csv": data_file.read_text(),
    }
