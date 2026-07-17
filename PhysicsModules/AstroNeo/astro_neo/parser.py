"""
Input-file parsing for Astro Neo: .ini -> validated flat parameter dict,
mirroring the NanoIndentation/EXAFS parser pattern (plain functions, no
import-time argparse side effects).

Unlike the original astro_neo_test copy, there is no [Paths] npaths/center/
fits section -- AstroNeo always fits exactly one hardcoded model (model.py),
so those keys never meant anything here; they were carried over unedited
from the EXAFS/XPS .ini schema this was copy-pasted from.
"""

import configparser

from PhysicsModules.AstroNeo.astro_neo.helper import bcolors, str_to_bool


def CheckKey(section, key_list):
    for key in key_list:
        if key not in section:
            raise KeyError(f"{key} is missing")


def optional_var(section, name_var, alt_var=None, type_var=int):
    if type_var == bool:
        if name_var in section:
            return str_to_bool(section[name_var])
        return alt_var
    if type_var is None:
        return section.get(name_var)
    if name_var in section:
        return type_var(section[name_var])
    return type_var(alt_var)


def read_input_file(input_file, verbose=False):
    """Read an Astro Neo .ini file into a dict of raw sections."""
    config = configparser.ConfigParser()
    read = config.read(input_file)
    if not read:
        raise FileNotFoundError(f"Input file not found: {input_file}")
    config = config._sections

    file_min = ["Inputs", "Populations", "Mutations", "Outputs"]
    CheckKey(config, file_min)

    CheckKey(config["Inputs"], ["data_dir", "data_file", "output_file"])
    CheckKey(config["Populations"], ["population", "num_gen"])
    CheckKey(config["Mutations"], ["mutated_options"])
    CheckKey(config["Outputs"], ["print_graph"])

    if verbose:
        print_input_file(config)
    return config


def validate_input_file(file_dict):
    """Flatten raw .ini sections into typed Astro Neo parameters."""
    inputs = file_dict["Inputs"]
    populations = file_dict["Populations"]
    mutations = file_dict["Mutations"]
    outputs = file_dict["Outputs"]

    return {
        # Inputs
        "data_dir": inputs["data_dir"],
        "data_file": inputs["data_file"],
        "bg_file": inputs.get("bg_file"),
        "rsp_file": inputs.get("rsp_file"),
        "acx2_path": inputs.get("acx2_path"),
        "output_file": inputs["output_file"],
        "xmin": optional_var(inputs, "xmin", 7.0, float),
        "xmax": optional_var(inputs, "xmax", 30.0, float),
        # Populations
        "nPops": int(populations["population"]),
        "nGen": int(populations["num_gen"]),
        "nBestSample": optional_var(populations, "best_sample", 0.3, float),
        "nLuckSample": optional_var(populations, "lucky_few", 0.2, float),
        # Mutations / solver
        "solOpt": optional_var(mutations, "solver_type", 2, int),
        "mutOpt": optional_var(mutations, "mutated_options", 4, int),
        "croOpt": optional_var(mutations, "crossover_options", 6, int),
        "mutChance": optional_var(mutations, "chance_of_mutation", 0.3, float),
        "F": optional_var(mutations, "mutf", 0.5, float),
        "CR": optional_var(mutations, "mutcr", 0.9, float),
        # Outputs
        "printGraph": optional_var(outputs, "print_graph", False, bool),
        "distributed": optional_var(outputs, "distributed", 1, int),
    }


def parse_input_file(input_file, verbose=False):
    """Convenience wrapper: .ini path -> validated flat parameter dict."""
    return validate_input_file(read_input_file(input_file, verbose=verbose))


def print_input_file(file_dict):
    for key, value in file_dict.items():
        print("[" + bcolors.BOLD + str(key) + bcolors.ENDC + "]")
        for inner_key, inner_value in value.items():
            print("---" + inner_key + ": " + str(inner_value))
