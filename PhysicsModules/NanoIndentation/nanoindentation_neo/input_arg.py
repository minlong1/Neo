import argparse
import sys

from PhysicsModules.NanoIndentation.nanoindentation_neo.nano_indent import NanoNeo
from PhysicsModules.NanoIndentation.nanoindentation_neo.parser import read_input_file


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", help="Submit input file to Nano Neo")
    parser.add_argument("-v", "--verbose", help="output verbosity", action="store_true")
    parser.add_argument("-s", "--show_input", help="show input file", action="store_true")

    args = parser.parse_args()
    if args.input is None:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.show_input:
        read_input_file(args.input, verbose=True)

    nano_neo = NanoNeo()
    nano_neo.nano_read(filepath=args.input)
    nano_neo.nano_setup()
    nano_neo.run()


if __name__ == "__main__":
    main()
