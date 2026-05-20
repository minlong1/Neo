import copy
import os
import subprocess
from pathlib import Path

from PhysicsModules.EXAFS.exafs_neo.utils import checkKey, check_if_exists

from typing import Dict, AnyStr, Union


def check_output_files(file: AnyStr):
    """
    Check if the output file for each file exists, else, assumes equality.

    Checks:
        <file>.csv
        <file>_data.csv
        <file>_generations.csv
    """
    check_if_exists(file)
    # self.file = file

    file_initial = open(file, "a+")
    file_initial.write(
        "Gen,TPS,FITTNESS,CO_Score,Mut_Score,CURRFIT,CURRIND,BESTFIT,BESTIND\n"
    )  # writing header
    file_initial.close()

    file_data = os.path.splitext(file)[0] + "_data.csv"
    check_if_exists(file_data)

    # Not using right now
    file_gen = os.path.splitext(file)[0] + "_generations.csv"
    check_if_exists(file_gen)


def sabcor_executable(
    base: Union[Path, str],
    sabcor_file: Union[Path, str],
    data_path: Union[Path, str],
    logger=None,
) -> Union[Path, str]:
    """
    Call the sabcor executable from the submodules

    The sabcor exectuable is as follows:

    sabcor <file.test> <inp>

    if input file is <file.test>, output is <file_sac.test>

    """

    # TODO: need to check if sabcor has been compiled
    if sabcor_file is not None:
        sabcor_exec = os.path.join(base, "contrib/sabcor/bin/sabcor")
        sabcor_full_file = os.path.join(base, sabcor_file)
        # call Sabcor
        command = [sabcor_exec, data_path, sabcor_full_file]
        p = subprocess.call(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        if p == 0:
            if logger is not None:
                logger.print("Sabcor Finished")
            return os.path.splitext(data_path)[0] + "_sac.csv"


class ExafsFileobj:
    sabcorFile: str

    def __init__(self):
        self.sabcor_toggle = None
        self.end = None
        self.base = None
        self._front = None
        self.nComp = 0
        self.data_file = None
        self.csv_series = False
        self.feff_file = None
        self.output_file = None

        self.data_path = None
        self.output_path = None
        self.log_path = None
        """

        Args:
            i (int, optional): the i-th file. Defaults to 0.
            firstPass (bool, optional): if is the first pass through the dataset. Defaults to False.
            path_optimize (bool, optional): if path optimize. Defaults to False.
        """

    def initialize(
        self,
        data_dict: Dict,
        i=0,
        firstPass: bool = False,
        path_optimize: bool = False,
        logger=None,
    ):
        """
        Initialize file paths for each of the file first

        @param dict data_dict:
        @param int i: the i-th file. Defaults to 0.
        @param bool firstPass: if is the first pass through the dataset. Defaults to False.
        @param bool path_optimize: if path optimization
        @param NeoLogger logger: logger
        @return:
        """
        # self.csv_series = csv_series
        self.base = os.getcwd()

        self.nComp = checkKey("nComp", data_dict, 1, logger=logger)
        self.data_file = checkKey("data_file", data_dict, logger=logger)
        self.feff_file = checkKey("feff_file", data_dict, logger=logger)
        self.output_file = checkKey("output_file", data_dict, logger=logger)
        self.sabcor_toggle = checkKey("sabcor_toggle", data_dict, logger=logger)
        self.sabcorFile = checkKey("sabcor_file", data_dict, False, logger=logger)

        if self.nComp > 1:
            self._front = []
            for i in range(self.nComp):
                self._front.append(os.path.join(self.base, self.feff_file[i]))
        else:
            self._front = os.path.join(self.base, self.feff_file)

        if self.csv_series:
            self.data_path = os.path.join(self.base, self.data_file[i])
            self.output_path = (
                os.path.splitext(os.path.join(self.base, self.output_file))[0]
                + "_"
                + str(i)
                + ".csv"
            )
            self.log_path = (
                os.path.splitext(copy.deepcopy(self.output_path))[0] + ".log"
            )

        else:
            self.data_path = os.path.join(self.base, self.data_file)
            self.output_path = os.path.join(self.base, self.output_file)
            self.log_path = (
                os.path.splitext(copy.deepcopy(self.output_path))[0] + ".log"
            )
        if path_optimize:
            self.output_path = (
                os.path.splitext(os.path.join(self.base, self.output_file))[0]
                + "_optimized.csv"
            )

        self.end = ".dat"
        check_output_files(self.output_path)
        if not firstPass:
            check_if_exists(self.log_path)

        if self.sabcor_toggle:
            sabcor_executable(self.base, self.sabcorFile, self.data_file, logger=logger)


if __name__ == "__main__":
    exafs_File = ExafsFileobj()

    data_dict = {
        "data_file": "path_files/Pu_C/pu3in_t030_sac.chi",
        "feff_file": ["path_files/Pu_C/feff/feff", "path_files/Pu_C/feff/feff2"],
        "output_file": "result/test/PuC.csv",
        "nComp": 2,
        "sabcor_toggle": False,
    }

    exafs_File.initialize(data_dict)

    print(exafs_File.output_file)
