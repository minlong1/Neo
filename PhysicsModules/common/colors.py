"""
Shared ANSI terminal color codes.

EXAFS (Bcolors, STRColors) and NanoIndentation (bcolors) each redefined the
same eight codes; this is the one canonical copy.
"""


class TermColors:
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"
