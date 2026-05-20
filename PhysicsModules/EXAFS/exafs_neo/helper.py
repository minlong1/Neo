import time
from typing import AnyStr

import numpy as np

from PhysicsModules.EXAFS.exafs_neo._version import __version__


class Bcolors:
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"


def time_call():
    return time.time()


def str_to_bool(s: AnyStr):
    s = s.lower()
    if s == "true":
        return True
    elif s == "false":
        return False
    else:
        raise ValueError  # evil ValueError that doesn't tell you what the wrong value was


def norm(val: float):
    return np.linalg.norm(val)


def banner():
    banner_str = r"""
            EXAFS_GA ver %s
 _______________________________________
|    _______  __    _    _____ ____     |
|   | ____\ \/ /   / \  |  ___/ ___|    |
|   |  _|  \  /   / _ \ | |_  \___ \    |
|   | |___ /  \  / ___ \|  _|  ___) |   |
|   |_____/_/\_\/_/   \_\_|   |____/    |
|_______________________________________|
    """ % __version__

    return banner_str
