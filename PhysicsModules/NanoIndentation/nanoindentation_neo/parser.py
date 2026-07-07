"""
Input-file parsing for Nano Neo.

Replaces the original module-level ini_parser/parser pair (which executed at
import time against global argparse state) with plain functions: .ini file
in, validated flat parameter dict out. The .ini sections and keys are
unchanged from nano_indent.
"""

import configparser

from PhysicsModules.NanoIndentation.nanoindentation_neo.helper import bcolors, str_to_bool


def CheckKey(section, key_list):
    for key in key_list:
        if key not in section:
            raise KeyError(f"{key} is missing")


def split_string(arr_str):
    """
    Read a comma-separated (optionally bracketed) list of floats
    """
    starter = [i for i, c in enumerate(arr_str) if c == "["]
    end = [i for i, c in enumerate(arr_str) if c == "]"]
    assert len(starter) == len(end), "Bracket setup not right."

    if starter:
        arr_str = arr_str[starter[0] + 1 : end[0]]
    return [float(i) for i in arr_str.split(",")]


def optional_var(section, name_var, alt_var=None, type_var=int):
    """
    Detections of optional variables exists within input files, and
        put in corresponding default inputs parameters.
    """
    if type_var == bool:
        if name_var in section:
            return str_to_bool(section[name_var])
        return alt_var
    if type_var is None:
        return section.get(name_var)
    if name_var in section:
        return type_var(section[name_var])
    return type_var(alt_var)


def optional_range(section, label):
    if label not in section:
        return []
    return split_string(section[label])


def read_input_file(input_file, verbose=False):
    """Read a Nano Neo .ini file into a dict of raw sections."""
    config = configparser.ConfigParser()
    read = config.read(input_file)
    if not read:
        raise FileNotFoundError(f"Input file not found: {input_file}")
    config = config._sections

    file_min = ["Inputs", "Populations", "Mutations", "Paths", "Outputs"]
    CheckKey(config, file_min)

    CheckKey(config["Inputs"], ["data_file", "output_file", "data_cutoff"])
    CheckKey(config["Populations"], ["population", "num_gen", "best_sample", "lucky_few"])
    CheckKey(
        config["Mutations"],
        ["chance_of_mutation", "original_chance_of_mutation", "mutated_options"],
    )
    CheckKey(config["Paths"], ["npaths", "fits"])
    CheckKey(config["Outputs"], ["print_graph", "num_output_paths"])

    if verbose:
        print_input_file(config)

    return config


def validate_input_file(file_dict):
    """Flatten raw .ini sections into typed Nano Neo parameters."""
    inputs = file_dict["Inputs"]
    populations = file_dict["Populations"]
    mutations = file_dict["Mutations"]
    paths = file_dict["Paths"]
    outputs = file_dict["Outputs"]

    return {
        # Inputs
        "data_file": inputs["data_file"],
        "output_file": inputs["output_file"],
        "data_cutoff": [float(x) for x in inputs["data_cutoff"].split(",")],
        # Populations
        "nPops": int(populations["population"]),
        "nGen": int(populations["num_gen"]),
        "nBestSample": int(populations["best_sample"]) / 100.0,
        "nLuckSample": int(populations["lucky_few"]) / 100.0,
        # Mutations (percent in the file, probability internally)
        "mutChance": int(mutations["chance_of_mutation"]) / 100.0,
        "mutOpt": int(mutations["mutated_options"]),
        # Paths
        "npaths": int(paths["npaths"]),
        "fits": paths["fits"].strip(),
        "A_range": optional_range(paths, "a_range"),
        "hf_range": optional_range(paths, "hf_range"),
        "m_range": optional_range(paths, "m_range"),
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
            print("---" + inner_key + ": " + inner_value)
