from PhysicsModules.EXAFS.exafs_neo._version import __version__
from PhysicsModules.common.cli import str_to_bool  # noqa: F401 (re-exported)
from PhysicsModules.common.colors import TermColors as Bcolors  # noqa: F401 (re-exported)


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
