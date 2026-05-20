import argparse
import os
import sys

from PhysicsModules.EXAFS.exafs_neo.exafs import ExafsNeo
from PhysicsModules.EXAFS.exafs_neo.ini_parser import validate_input_file
from PhysicsModules.EXAFS.exafs_neo.parser import InputParamsParser

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help="Submit input file to EXAFS")
parser.add_argument("-v", "--verbose", help="output verbosity", action="store_true")
parser.add_argument("-s", "--show_input", help="show input file", action="store_true")
parser.add_argument("-t", help="Timeing mode", action="store_true")
parser.add_argument("-d", help="Debug mode", action="store_true")

args = parser.parse_args()
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

if args.input is not None:
    file_path = os.path.join(os.getcwd(), args.input)
    input_params = InputParamsParser()
    input_params.read_input_file(file_path, verbose=args.show_input)
    input_params.input_dict = validate_input_file(input_params.input_dict)

    input_pars = input_params.export_input_dict()

else:
    print("No input file is given")

debug_mode = args.d
timeing_mode = args.t


def main():
    exafs_neo = ExafsNeo()
    exafs_neo.exafs_read(input_parameters=input_pars)
    exafs_neo.exafs_setup()
    exafs_neo.run()


if __name__ == "__main__":
    main()
