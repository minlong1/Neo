"""
Shared CLI/input-file helpers used by every physics module: boolean parsing
for .ini values and a minimal stdout(+file) logger. Each module's own
NeoLogger/NanoLogger differed only in logger name and default level, so
BaseLogger takes those as parameters instead.
"""

import logging
import sys


def str_to_bool(s: str) -> bool:
    normalized = str(s).strip().lower()
    if normalized == "true":
        return True
    if normalized == "false":
        return False
    raise ValueError(f"Expected 'True' or 'False', got {s!r}")


class BaseLogger:
    def __init__(self, name: str = "", level: int = logging.INFO):
        logging.basicConfig(level=level)
        self.logger = logging.getLogger(name)
        if self.logger.hasHandlers():
            self.logger.handlers.clear()
        self.log_path = None
        self.logging_level = level

    def initialize_logging(self, log_path=None, log_format: str = "%(message)s") -> None:
        self.log_path = log_path
        formatter = logging.Formatter(log_format)

        if log_path is not None:
            file_handler = logging.FileHandler(log_path, mode="a+", encoding="utf-8")
            file_handler.setFormatter(formatter)
            file_handler.setLevel(self.logging_level)
            self.logger.addHandler(file_handler)

        stdout_handler = logging.StreamHandler(stream=sys.stdout)
        stdout_handler.setLevel(self.logging_level)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(stdout_handler)
        self.logger.setLevel(self.logging_level)

    def set_loglevel(self, loglevel: int) -> None:
        self.logging_level = loglevel
        self.logger.setLevel(loglevel)

    def print(self, message: str) -> None:
        self.logger.log(self.logging_level, message)

    def __call__(self, message: str) -> None:
        self.print(message)
