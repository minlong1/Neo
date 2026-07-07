# Parser for Inputs files
"""
Renamed to parser2 because it was overshadowing the stdlib parser when importing it non-relative for dev stuff - Evan 7/7/23
"""
import configparser
from typing import Dict, List

from PhysicsModules.XPS.xps_neo.helper import bcolors

def CheckKey(key_dict: Dict, key_list: List) -> None:
    for i in range(len(key_list)):
        try:
            key_dict[key_list[i]]
        except KeyError:
            raise KeyError(str(key_list[i]) + ' is missing')
            # break

def CheckOptionalKey(key_dict, optional_key_list):
    optional_key = []
    for i in range(len(optional_key_list)):
        try:
            key_dict[optional_key_list[i]]
        except KeyError:
            optional_key.append(optional_key_list[i])

    return optional_key
# config the parser to the right order:
def read_input_file(input_file,verbose=False):
    config = configparser.ConfigParser()
    config.read(input_file)
    config = config._sections

    # read into each dict
    file_min = ['Inputs','Populations','Mutations','Paths','Outputs']
    CheckKey(config,file_min)

    Inputs_dict = config['Inputs']
    Populations_dict = config['Populations']
    Mutations_dict = config['Mutations']
    Paths_dict = config['Paths']
    Outputs_dict = config['Outputs']

    # Checking for minimum inputs
    #Removed 'data_cutoff' -evan
    input_min = ['data_file','output_file']
    CheckKey(Inputs_dict,input_min)

    population_min = ['population','num_gen','best_sample','lucky_few']
    CheckKey(Populations_dict,population_min)

    mutation_min = ['chance_of_mutation','original_chance_of_mutation','mutated_options']
    CheckKey(Mutations_dict,mutation_min)

    path_min = ['npeaks']
    path_optional = ['path_optimize','steady_state','corr_paths']
    #print(Paths_dict.keys())
    CheckKey(Paths_dict,path_min)
    path_missing = CheckOptionalKey(Paths_dict,path_optional)


    # larch_min = ['kmin','kmax','kweight','deltak','rbkg','bkgkw','bkgkmax']
    # CheckKey(Larch_dict,larch_min)

    output_min =['print_graph','num_output_paths']
    output_optional = ['steady_state']
    CheckKey(Outputs_dict,output_min)
    output_missing = CheckOptionalKey(Outputs_dict,output_min)

    # Adjust values

    # Pack all of them into a single dicts
    file_dict = {}
    file_dict['Inputs'] = Inputs_dict
    file_dict['Populations'] = Populations_dict
    file_dict['Mutations'] = Mutations_dict
    file_dict['Paths'] = Paths_dict
    file_dict['Outputs'] = Outputs_dict

    if verbose == True:
        print_input_file(file_dict)

    return file_dict

def print_input_file(file_dict):
    for key,value in file_dict.items():
        print("[" +bcolors.BOLD + str(key)+ bcolors.ENDC +"]")
        for inner_key,inner_value in value.items():
            print('---'  +  inner_key + ": " +inner_value)


## Need to run some sample:
# if __name__ == '__main__':
#
#     checkKey()
