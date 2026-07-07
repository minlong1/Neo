import logging

from PhysicsModules.common.cli import BaseLogger, str_to_bool  # noqa: F401 (re-exported)
from PhysicsModules.common.colors import TermColors as bcolors  # noqa: F401 (re-exported)

__version__ = "0.1.0"


class NanoLogger(BaseLogger):
    def __init__(self):
        super().__init__(name="nano_neo", level=logging.INFO)


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
