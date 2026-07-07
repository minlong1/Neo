import logging
import sys

__version__ = "0.1.0"


class bcolors:
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"


def str_to_bool(s):
    if s == "True":
        return True
    elif s == "False":
        return False
    else:
        raise ValueError(f"Expected 'True' or 'False', got {s!r}")


class NanoLogger:
    """Message logger writing to stdout and, optionally, a log file."""

    def __init__(self):
        self.logger = logging.getLogger("nano_neo")
        self.logger.handlers = []
        self.logger.setLevel(logging.INFO)
        self.log_path = None

    def initialize_logging(self, log_path=None, log_format="%(message)s"):
        formatter = logging.Formatter(log_format)
        self.log_path = log_path

        if log_path is not None:
            file_handler = logging.FileHandler(log_path, mode="a+", encoding="utf-8")
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)

        stdout_handler = logging.StreamHandler(stream=sys.stdout)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)

    def print(self, message: str):
        self.logger.info(message)

    def __call__(self, message):
        self.logger.info(message)


def banner():
    """
    https://patorjk.com/software/taag/#p=display&h=1&v=1&f=Big%20Chief&t=Nano%20Neo
    """
    banner_str = (
        """
                    Nano_Neo ver %s
_____________________________________________________________
    _     _                             _     _
    /|   /                              /|   /
---/-| -/------__-----__-----__--------/-| -/------__-----__-
  /  | /     /   )  /   )  /   )      /  | /     /___)  /   )
_/___|/_____(___(__/___/__(___/______/___|/_____(___ __(___/_

    """
        % __version__
    )

    return banner_str
