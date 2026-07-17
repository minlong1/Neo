import logging

from PhysicsModules.common.cli import BaseLogger, str_to_bool  # noqa: F401 (re-exported)
from PhysicsModules.common.colors import TermColors as bcolors  # noqa: F401 (re-exported)

__version__ = "0.1.0"


class AstroLogger(BaseLogger):
    def __init__(self):
        super().__init__(name="astro_neo", level=logging.INFO)


def banner():
    banner_str = (
        """
            Astro-Neo ver %s
     _        _               _   _
    / \\   ___| |_ _ __ ___   | \\ | | ___  ___
   / _ \\ / __| __| '__/ _ \\  |  \\| |/ _ \\/ _ \\
  / ___ \\\\__ \\ |_| | | (_) | | |\\  |  __/ (_) |
 /_/   \\_\\___/\\__|_|  \\___/  |_| \\_|\\___|\\___/
"""
        % __version__
    )
    return banner_str
