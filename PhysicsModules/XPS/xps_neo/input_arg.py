from PhysicsModules.XPS.xps_neo import __version__
from PhysicsModules.XPS.xps_neo.parser2 import read_input_file
import argparse, os, sys

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help="Submit input file")
parser.add_argument("--version", action="version", version=f"xps_neo {__version__}")
parser.add_argument("-v", "--verbose", help="output verbosity", action="store_true")
parser.add_argument("-s", "--show_input", help="show input file", action="store_true")
parser.add_argument("-t", help="Timeing mode", action="store_true")
parser.add_argument("-d", help="Debug mode", action="store_true")
parser.add_argument(
    "--seed",
    type=int,
    default=None,
    help="Seed the random number generators for a reproducible run",
)
parser.add_argument(
    "--workers",
    type=int,
    default=1,
    help="Parallel processes for population evaluation (default 1 = "
    "serial; results are identical either way)",
)

# Import-time defaults: modules that star-import this one (import_lib) may
# consult these before parse_args() has run. They are updated by parse_args.
debug_mode = False
timeing_mode = False


def parse_args(argv=None):
    """Parse CLI arguments and apply the RNG seed. No work at import time."""
    global debug_mode, timeing_mode
    if argv is None:
        argv = sys.argv[1:]
    if len(argv) == 0:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args(argv)

    if args.seed is not None:
        import random
        import numpy as np

        random.seed(args.seed)
        np.random.seed(args.seed)

    debug_mode = args.d
    timeing_mode = args.t
    return args


def load_file_dict(args):
    """Read the INI named by -i/--input into a section dict."""
    if args.input is None:
        parser.print_help(sys.stderr)
        sys.exit(1)
    file_path = os.path.join(os.getcwd(), args.input)
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"input file not found: {file_path}")
    file_dict = read_input_file(file_path)
    if args.show_input:
        read_input_file(file_path, verbose=True)
    return file_dict
