from PhysicsModules.EXAFS.exafs_neo.helper import str_to_bool
from typing import List, Dict

"""
Author: Andy Lau
Last Updated: 1/19/2021

Changes:

2/8/2021: Andy:
    - Function definition for arr to be in 2D

1/20/2021: Andy:
    - Function definition for optional parameters.
    - Changes to allows for multiple inputs folder.

"""


def split_path_arr(arr_str: str, num_comps: int) -> List:
    """
    Read the path list

    @param str arr_str: str of array for the path list
    @param int num_comps: number of compounds
    @return:
    """
    assert num_comps > 0
    starter = []
    end = []
    k = 0
    split_str = []
    for i in arr_str:
        if i == "[":
            starter.append(k)
        elif i == "]":
            end.append(k)
        k = k + 1

    assert len(starter) == len(end), "Bracket setup not matching."
    if num_comps > 1:
        assert num_comps == len(starter), "Number of compounds not matched."
        assert num_comps == len(end), "Number of compounds not matched."

    # check if both are zeros, therefore the array is one 1 dimensions
    if len(starter) == 0 and len(end) == 0:
        split_str = list(arr_str.split(","))
    else:
        for i in range(len(starter)):
            split_str.append(arr_str[starter[i] + 1 : end[i]].split(","))

    return split_str


def optional_var(
    input_dict: Dict, name_var: str, alt_var=None, type_var=int, output_var: bool = True
):
    """
    Detections of optional variables exists within input files, and put in corresponding default inputs parameters.
    boolean needs special attentions
    """
    if type_var == bool:
        if name_var in input_dict:
            return_var = str_to_bool(input_dict[name_var])
        else:
            return_var = alt_var
    elif type_var is None:
        if name_var in input_dict:
            return_var = input_dict[name_var]
        else:
            return_var = None
    else:
        if name_var in input_dict:
            return_var = type_var(input_dict[name_var])
        else:
            return_var = type_var(alt_var)

    # return return_var
    input_dict[name_var] = return_var
    if output_var:
        return return_var


def validate_input_file(file_dict: Dict) -> Dict:
    Inputs_dict = file_dict["Inputs"]
    Populations_dict = file_dict["Populations"]
    Solver_dict = file_dict["Solver"]
    Mutations_dict = file_dict["Mutations"]
    Paths_dict = file_dict["Paths"]
    Larch_dict = file_dict["Larch_Paths"]
    Outputs_dict = file_dict["Outputs"]

    # Inputs
    num_compounds = optional_var(Inputs_dict, "num_compounds", 1, int)
    csv_file = Inputs_dict["csv_file"]
    output_file = Inputs_dict["output_file"]
    pathrange_file = optional_var(Inputs_dict, "pathrange_file", None, None)
    sabcor_file = optional_var(Inputs_dict, "sabcor_file", None, None)

    # Compounds
    if num_compounds > 1:
        try:
            feff_file = list(Inputs_dict["feff_file"].split(","))
        except:
            print("Feff folder is not correct")
    else:
        feff_file = Inputs_dict["feff_file"]

    try:
        csv_series = str_to_bool(Inputs_dict["csv_series"])
        if csv_series:
            csv_file = list(Inputs_dict["csv_file"].split(","))
    except KeyError:
        csv_series = False

    # population
    size_population = int(Populations_dict["population"])
    number_of_generation = int(Populations_dict["num_gen"])
    Populations_dict["best_sample"] = (
        float(Populations_dict["best_sample"]) / size_population
    )
    Populations_dict["lucky_few"] = (
        float(Populations_dict["lucky_few"]) / size_population
    )
    # Solver:
    solver_options = optional_var(Solver_dict, "solver_options", 0, int)

    # Mutations
    Mutations_dict["chance_of_mutation"] = 0.01 * float(
        Mutations_dict["chance_of_mutation"]
    )
    Mutations_dict["original_chance_of_mutation"] = 0.01 * float(
        (Mutations_dict["original_chance_of_mutation"])
    )
    Mutations_dict["chance_of_mutation_e0"] = 0.01 * float(
        (Mutations_dict["chance_of_mutation_e0"])
    )
    selection_options = optional_var(Mutations_dict, "selection_options", 0, int)
    mutation_options = optional_var(Mutations_dict, "mutated_options", 0, int)
    crossover_options = optional_var(Mutations_dict, "crossover_options", 0, int)
    # mutated_options = int(Mutations_dict['mutated_options'])

    # Paths
    if num_compounds > 1:
        individual_path = True
    else:
        individual_path = str_to_bool(Paths_dict["individual_path"])
        pathrange = int(Paths_dict["path_range"])

    try:
        optimize_only = str_to_bool(Paths_dict["optimize_only"])
    except KeyError:
        optimize_only = False

    try:
        path_optimize = str_to_bool(Paths_dict["path_optimize"])
    except KeyError:
        path_optimize = False

    try:
        Paths_dict["path_optimize_percent"] = float(Paths_dict["path_optimize_percent"])
    except KeyError:
        Paths_dict["path_optimize_percent"] = 0.01

    Paths_dict["path_list"] = split_path_arr(Paths_dict["path_list"], num_compounds)
    Paths_dict["path_optimize"] = optional_var(Paths_dict, "path_optimize", False, bool)
    Paths_dict["path_optimize_percent"] = path_optimize_percent = optional_var(
        Paths_dict, "path_optimize_percent", 0.01, float
    )
    # Larch Paths
    kmin = float(Larch_dict["kmin"])
    kmax = float(Larch_dict["kmax"])
    kweight = float(Larch_dict["kweight"])
    deltak = float(Larch_dict["deltak"])
    rbkg = float(Larch_dict["rbkg"])
    bkgkw = float(Larch_dict["bkgkw"])
    bkgkmax = float(Larch_dict["bkgkmax"])

    # Output
    printgraph = str_to_bool(Outputs_dict["print_graph"])
    num_output_paths = str_to_bool(Outputs_dict["num_output_paths"])
    steady_state = optional_var(Outputs_dict, "steady_state_exit", False, bool)

    # Map it back into single dictionary
    temp_dict = {
        "Inputs": Inputs_dict,
        "Populations": Populations_dict,
        "Solver": Solver_dict,
        "Mutations": Mutations_dict,
        "Paths": Paths_dict,
        "Larch_Paths": Larch_dict,
        "Outputs": Outputs_dict,
    }
    # Package into a single dictionary
    return temp_dict
