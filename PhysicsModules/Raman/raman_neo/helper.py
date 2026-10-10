import logging

from PhysicsModules.common.cli import BaseLogger, str_to_bool  # noqa: F401 (re-exported)
from PhysicsModules.common.colors import TermColors as bcolors  # noqa: F401 (re-exported)

# Historical raman_neo name for the color codes
Bcolors = bcolors

__version__ = "0.1.0"


class RamanLogger(BaseLogger):
    def __init__(self):
        super().__init__(name="raman_neo", level=logging.INFO)


def banner():
    banner_str = (
        r"""
            Raman Neo ver %s
 _________________________________________
|   ____      _    __  __    _    _   _    |
|  |  _ \    / \  |  \/  |  / \  | \ | |   |
|  | |_) |  / _ \ | |\/| | / _ \ |  \| |   |
|  |  _ <  / ___ \| |  | |/ ___ \| |\  |   |
|  |_| \_\/_/   \_\_|  |_/_/   \_\_| \_|   |
|__________________________________________|
"""
        % __version__
    )
    return banner_str
