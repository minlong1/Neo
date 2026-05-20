"""
Authors    Matthew Adas, Miu Lun(Andy) Lau*, Jeffrey Terry, Min Long
Email      madas@hawk.iit.edu, andylau@u.boisestate.edu
Version    0.2.0
Date       July 4, 2021

Constains misc functions for GUI
"""

import platform
import re
import sys
import os


def get_platform() -> str:
    os_name = platform.system()
    return os_name


def atoi(text) -> None:
    return int(text) if text.isdigit() else text


def natural_keys(text):
    """
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    """
    return [atoi(c) for c in re.split(r"(\d+)", text)]


def check_sabcor_folder(ostype: str):
    # check sabcor folder if it initializes
    default_sabcor = False
    if ostype == "Windows":
        print("Window doesn't have sabcor capabilites, continue without sabcor")
    else:
        if os.path.exists("./contrib/sabcor"):
            if not os.listdir("./contrib/sabcor"):
                print('Please Initialize "Sabcor" as submodules')
                exit()
            else:
                if os.path.exists("../contrib/sabcor/sabcor"):
                    print(
                        "Please make the sabcor executable first in './contrib/sabcor' directory"
                    )
                else:
                    sys.path.insert(1, "./contrib/sabcor")

                    default_sabcor = True
        else:
            print("Sabcor is not found...")

    return default_sabcor


if __name__ == "__main__":
    os_name = get_platform()
    print(os_name)
