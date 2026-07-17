import argparse
import sys

from PhysicsModules.AstroNeo.astro_neo.astro_neo import AstroNeo
from PhysicsModules.AstroNeo.astro_neo.parser import read_input_file


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", help="Submit input file to Astro Neo")
    parser.add_argument("-v", "--verbose", help="output verbosity", action="store_true")
    parser.add_argument("-s", "--show_input", help="show input file", action="store_true")

    args = parser.parse_args()
    if args.input is None:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.show_input:
        read_input_file(args.input, verbose=True)

    astro_neo = AstroNeo()
    astro_neo.astro_read(filepath=args.input)
    astro_neo.astro_setup()
    astro_neo.run()


if __name__ == "__main__":
    main()
