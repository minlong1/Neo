"""
Input-file parsing for Raman Neo.

Replaces raman_neo's InputParamsParser/ini_parser pair with plain functions
(the NanoIndentation/AstroNeo pattern): .ini file in, validated flat
parameter dict out. The .ini sections and keys are unchanged from the
standalone raman_neo package:

    [Inputs]      data_file, output_file, band_file,
                  optional num_compounds, bandrange_file
    [Populations] population, num_gen, best_sample, lucky_few
                  (best_sample/lucky_few are counts, not percents)
    [Solver]      optional solver_options (0 GA, 1 GA+Rechenberg, 2 DE, or a
                  Solvers registry name such as PSO / CMA_ES), mutf, mutcr
    [Mutations]   chance_of_mutation (percent), optional mutated_options,
                  selection_options, crossover_options. original_chance_of_mutation
                  and chance_of_mutation_e0 are still accepted but unused: the
                  shift is an ordinary gene, mutated by the selected operator
    [Bands]       band_list ("1,2,3" or "[1,2][3,4]" for multiple compounds),
                  optional band_optimize, band_optimize_percent
    [Spectrum]    xmin, xmax, weight, deltax, optional baseline_order,
                  normalize, shift_range ("low, high, step")
    [Outputs]     print_graph, optional steady_state_exit
"""

import configparser

from PhysicsModules.Raman.raman_neo.helper import bcolors, str_to_bool

SECTIONS = ["Inputs", "Populations", "Solver", "Mutations", "Bands", "Spectrum", "Outputs"]


def CheckKey(section, key_list):
    for key in key_list:
        if key not in section:
            raise KeyError(f"{key} is missing")


def split_band_arr(arr_str: str, num_comps: int = 1):
    """Parse ``band_list``: ``"1,2,3"`` -> ``[1, 2, 3]``;
    ``"[1,2][3,4]"`` -> ``[[1, 2], [3, 4]]`` (one list per compound)."""
    if num_comps < 1:
        raise ValueError("num_compounds must be >= 1")
    starter = [i for i, c in enumerate(arr_str) if c == "["]
    end = [i for i, c in enumerate(arr_str) if c == "]"]
    if len(starter) != len(end):
        raise ValueError(f"Bracket setup not matching in band_list {arr_str!r}")
    if num_comps > 1 and len(starter) != num_comps:
        raise ValueError(
            f"band_list {arr_str!r} has {len(starter)} groups but num_compounds = {num_comps}"
        )

    def ints(chunk):
        return [int(float(v)) for v in chunk.split(",") if v.strip()]

    if not starter:
        return ints(arr_str)
    groups = [ints(arr_str[s + 1 : e]) for s, e in zip(starter, end)]
    return groups if num_comps > 1 else groups[0]


def split_floats(arr_str):
    """Comma-separated (optionally bracketed) floats."""
    arr_str = arr_str.strip().lstrip("[").rstrip("]")
    return [float(v) for v in arr_str.split(",") if v.strip()]


def optional_var(section, name_var, alt_var=None, type_var=int):
    """Typed optional key with a default (same helper as NanoIndentation)."""
    if type_var == bool:
        if name_var in section:
            return str_to_bool(section[name_var])
        return alt_var
    if type_var is None:
        return section.get(name_var, alt_var)
    if name_var in section:
        return type_var(section[name_var])
    return alt_var if alt_var is None else type_var(alt_var)


def solver_option(value):
    """``solver_options`` is a historical numeric ID (0/1/2) or a Solvers
    registry name (e.g. ``PSO``, ``CMA_ES``)."""
    value = str(value).strip()
    try:
        return int(value)
    except ValueError:
        return value


def read_input_file(input_file, verbose=False):
    """Read a Raman Neo .ini file into a dict of raw sections."""
    config = configparser.ConfigParser()
    read = config.read(input_file)
    if not read:
        raise FileNotFoundError(f"Input file not found: {input_file}")
    config = config._sections

    CheckKey(config, SECTIONS)
    CheckKey(config["Inputs"], ["data_file", "output_file", "band_file"])
    CheckKey(config["Populations"], ["population", "num_gen", "best_sample", "lucky_few"])
    CheckKey(config["Mutations"], ["chance_of_mutation"])
    CheckKey(config["Bands"], ["band_list"])
    CheckKey(config["Spectrum"], ["xmin", "xmax", "weight", "deltax"])
    CheckKey(config["Outputs"], ["print_graph"])

    if verbose:
        print_input_file(config)

    return config


def validate_input_file(file_dict):
    """Flatten raw .ini sections into typed Raman Neo parameters."""
    inputs = file_dict["Inputs"]
    populations = file_dict["Populations"]
    solver = file_dict["Solver"]
    mutations = file_dict["Mutations"]
    bands = file_dict["Bands"]
    spectrum = file_dict["Spectrum"]
    outputs = file_dict["Outputs"]

    num_compounds = optional_var(inputs, "num_compounds", 1, int)
    band_file = inputs["band_file"].strip()
    if num_compounds > 1:
        band_file = [f.strip() for f in band_file.split(",")]

    n_pops = int(populations["population"])
    shift_range = spectrum.get("shift_range")

    return {
        # Inputs
        "data_file": inputs["data_file"].strip(),
        "output_file": inputs["output_file"].strip(),
        "band_file": band_file,
        "bandrange_file": optional_var(inputs, "bandrange_file", None, None),
        "num_compounds": num_compounds,
        # Populations (counts in the file, fractions internally)
        "nPops": n_pops,
        "nGen": int(populations["num_gen"]),
        "nBestSample": float(populations["best_sample"]) / n_pops,
        "nLuckSample": float(populations["lucky_few"]) / n_pops,
        # Solver
        "solOpt": solver_option(solver.get("solver_options", 0)),
        "F": optional_var(solver, "mutf", 0.5, float),
        "CR": optional_var(solver, "mutcr", 0.9, float),
        # Mutations (percent in the file, probability internally)
        "mutChance": float(mutations["chance_of_mutation"]) / 100.0,
        "mutOpt": optional_var(mutations, "mutated_options", 1, int),
        "selOpt": optional_var(mutations, "selection_options", 0, int),
        "croOpt": optional_var(mutations, "crossover_options", 0, int),
        # Bands
        "band_list": split_band_arr(bands["band_list"], num_compounds),
        "bandOptimize": optional_var(bands, "band_optimize", False, bool),
        "bandOptimizePercent": optional_var(bands, "band_optimize_percent", 0.01, float),
        # Spectrum
        "xmin": float(spectrum["xmin"]),
        "xmax": float(spectrum["xmax"]),
        "weight": float(spectrum["weight"]),
        "deltax": float(spectrum["deltax"]),
        "baseline_order": optional_var(spectrum, "baseline_order", 0, int),
        "normalize": optional_var(spectrum, "normalize", True, bool),
        "shift_range": split_floats(shift_range) if shift_range else [],
        # Outputs
        "printGraph": str_to_bool(outputs["print_graph"]),
        "steadyState": optional_var(outputs, "steady_state_exit", False, bool),
    }


def parse_input_file(input_file, verbose=False):
    """Convenience wrapper: .ini path -> validated flat parameter dict."""
    return validate_input_file(read_input_file(input_file, verbose=verbose))


def print_input_file(file_dict):
    for key, value in file_dict.items():
        print("[" + bcolors.BOLD + str(key) + bcolors.ENDC + "]")
        for inner_key, inner_value in value.items():
            print("---" + inner_key + ": " + str(inner_value))
