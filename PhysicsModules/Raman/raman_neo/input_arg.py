import argparse
import sys

from PhysicsModules.Raman.raman_neo.parser import read_input_file
from PhysicsModules.Raman.raman_neo.raman import RamanNeo


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", help="Submit input file to Raman Neo")
    parser.add_argument("-v", "--verbose", help="output verbosity", action="store_true")
    parser.add_argument("-s", "--show_input", help="show input file", action="store_true")

    args = parser.parse_args()
    if args.input is None:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.show_input:
        read_input_file(args.input, verbose=True)

    raman_neo = RamanNeo()
    raman_neo.raman_read(filepath=args.input.strip())
    raman_neo.raman_setup()
    raman_neo.run()


if __name__ == "__main__":
    main()
