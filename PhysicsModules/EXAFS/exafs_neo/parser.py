# Parser for Inputs files
import configparser
from dataclasses import dataclass, field

from PhysicsModules.EXAFS.exafs_neo.helper import Bcolors


def CheckKey(dict_val, key_list):
    for i in range(len(key_list)):
        try:
            dict_val[key_list[i]]
        except KeyError:
            raise KeyError(str(key_list[i]) + " is missing")
            # break


def print_input_file(file_dict):
    for key, value in file_dict.items():
        print("[" + Bcolors.BOLD + str(key) + Bcolors.ENDC + "]")
        for inner_key, inner_value in value.items():
            print("---" + inner_key + ": " + inner_value)


def check_optional_key(og_dict, optional_key_list):
    optional_key = []
    for i in range(len(optional_key_list)):
        try:
            og_dict[optional_key_list[i]]
        except KeyError:
            optional_key.append(optional_key_list[i])

    return optional_key


@dataclass
class InputParamsParser:
    input_dict: dict = field(default_factory=dict)

    def read_input_file(self, input_file, verbose=False):
        config_parser = configparser.ConfigParser()
        config_parser.read(input_file)
        config = config_parser._sections
        # read into each dict
        file_min = [
            "Inputs",
            "Populations",
            "Solver",
            "Mutations",
            "Paths",
            "Larch_Paths",
            "Outputs",
        ]
        CheckKey(config, file_min)

        Inputs_dict = config["Inputs"]
        Populations_dict = config["Populations"]
        Solution_dict = config["Solver"]
        Mutations_dict = config["Mutations"]
        Paths_dict = config["Paths"]
        Larch_dict = config["Larch_Paths"]
        Outputs_dict = config["Outputs"]

        # Checking for minimum inputs
        input_min = ["csv_file", "output_file", "feff_file"]
        input_optional = ["num_compounds", "pathrange_file", "sabcor_file"]
        CheckKey(Inputs_dict, input_min)
        input_missing = check_optional_key(Inputs_dict, input_optional)

        population_min = ["population", "num_gen", "best_sample", "lucky_few"]
        CheckKey(Populations_dict, population_min)

        solver_optional = ["solver_options"]
        check_optional_key(Solution_dict, solver_optional)

        mutation_min = [
            "chance_of_mutation",
            "original_chance_of_mutation",
            "chance_of_mutation_e0",
        ]
        mutation_optional = [
            "mutated_options",
            "selection_options",
            "crossover_options",
        ]
        CheckKey(Mutations_dict, mutation_min)
        # mut_optional = CheckOptionalKey(Mutations_dict,mutation_optional)

        path_min = ["path_range", "path_list", "individual_path"]
        path_optional = ["path_optimize", "optimize_percent", "optimize_only"]
        CheckKey(Paths_dict, path_min)
        path_missing = check_optional_key(Paths_dict, path_optional)

        larch_min = ["kmin", "kmax", "kweight", "deltak", "rbkg", "bkgkw", "bkgkmax"]
        CheckKey(Larch_dict, larch_min)

        output_min = ["print_graph", "num_output_paths"]
        output_optional = ["steady_state_exit"]
        CheckKey(Outputs_dict, output_min)
        output_missing = check_optional_key(Outputs_dict, output_optional)
        # Adjust values

        # Pack all of them into a single dicts
        self.input_dict["Inputs"] = Inputs_dict
        self.input_dict["Populations"] = Populations_dict
        self.input_dict["Solver"] = Solution_dict
        self.input_dict["Mutations"] = Mutations_dict
        self.input_dict["Paths"] = Paths_dict
        self.input_dict["Larch_Paths"] = Larch_dict
        self.input_dict["Outputs"] = Outputs_dict

        if verbose:
            print_input_file(self.input_dict)

    def verbose(self):
        print_input_file(self.input_dict)

    def export_input_dict(self):
        """
        Convert the file dictionary from multiple dimensions into a singple dictionary
        """
        temp_dict = {
            # Input
            "num_compounds": self.input_dict["Inputs"]["num_compounds"],
            "data_file": self.input_dict["Inputs"]["csv_file"],
            "output_file": self.input_dict["Inputs"]["output_file"],
            "pathrange_file": self.input_dict["Inputs"]["pathrange_file"],
            # 'log_file': self.input_dict['Inputs']['log_file'],
            "feff_file": self.input_dict["Inputs"]["feff_file"],
            "sabcor_file": self.input_dict["Inputs"]["sabcor_file"],
            # Population
            "nPops": int(self.input_dict["Populations"]["population"]),
            "nGen": int(self.input_dict["Populations"]["num_gen"]),
            # Solver
            "solOpt": int(self.input_dict["Solver"]["solver_options"]),
            "selOpt": int(self.input_dict["Mutations"]["selection_options"]),
            "nBestSample": self.input_dict["Populations"]["best_sample"],
            "nLuckySample": int(self.input_dict["Populations"]["lucky_few"]),
            "steadyState": self.input_dict["Outputs"]["steady_state_exit"],
            "printGraph": self.input_dict["Outputs"]["print_graph"],
            "mut_options": int(self.input_dict["Mutations"]["mutated_options"]),
            "mutChance": self.input_dict["Mutations"]["chance_of_mutation"],
            "mutChanceE0": self.input_dict["Mutations"]["chance_of_mutation_e0"],
            "croOpt": int(self.input_dict["Mutations"]["crossover_options"]),
            # Larch
            "kmin": float(self.input_dict["Larch_Paths"]["kmin"]),
            "kmax": float(self.input_dict["Larch_Paths"]["kmax"]),
            "kweight": float(self.input_dict["Larch_Paths"]["kweight"]),
            "deltak": float(self.input_dict["Larch_Paths"]["deltak"]),
            "rbkg": float(self.input_dict["Larch_Paths"]["rbkg"]),
            "bkgkw": float(self.input_dict["Larch_Paths"]["bkgkw"]),
            "bkgkmax": float(self.input_dict["Larch_Paths"]["bkgkmax"]),
            # Paths
            "pathrange": self.input_dict["Paths"]["path_list"],
            "individualOptions": self.input_dict["Paths"]["individual_path"],
            "pathOptimize": self.input_dict["Paths"]["path_optimize"],
            "pathOptimize_percent": self.input_dict["Paths"]["path_optimize_percent"],
        }

        return temp_dict


## Need to run some sample:
# if __name__ == '__main__':
#
#     checkKey()
