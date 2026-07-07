import time

from PhysicsModules.common.cli import str_to_bool  # noqa: F401 (re-exported)
from PhysicsModules.common.colors import TermColors as bcolors  # noqa: F401 (re-exported)

from PhysicsModules.XPS.xps_neo.__init__ import __version__


def timecall():
    return time.time()


def str_to_list(s):
    return [float(i) for i in list(s.split(","))]


def banner():
    """
    https://patorjk.com/software/taag/#p=display&h=2&v=0&f=Univers&t=XPS%20NEO
    """
    banner_str = (
        '''
                                XPS Neo ver %s
__________________________________________________________________________________

8b        d8  88888888ba    ad88888ba      888b      88  88888888888  ,ad8888ba,
 Y8,    ,8P   88      "8b  d8"     "8b     8888b     88  88          d8"'    `"8b
  `8b  d8'    88      ,8P  Y8,             88 `8b    88  88         d8'        `8b
    Y88P      88aaaaaa8P'  `Y8aaaaa,       88  `8b   88  88aaaaa    88          88
    d88b      88""""""'      `"""""8b,     88   `8b  88  88"""""    88          88
  ,8P  Y8,    88                   `8b     88    `8b 88  88         Y8,        ,8P
 d8'    `8b   88           Y8a     a8P     88     `8888  88          Y8a.    .a8P
8P        Y8  88            "Y88888P"      88      `888  88888888888  `"Y8888Y"'
__________________________________________________________________________________

    '''
        % __version__
    )

    return banner_str
