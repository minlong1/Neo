from pathlib import Path
from typing import Dict

from attrs import define, field

from PhysicsModules.EXAFS.exafs_neo.utils import checkKey


@define
class NeoFilePars:
    # Default base to the EXAFS package directory so relative data paths in the
    # test suite (e.g. "tests/cu_test_files/...") resolve correctly even when
    # tests are discovered from the repository root or another cwd.
    base: Path = Path(__file__).resolve().parents[1]
    data_file: str = ""
    output_file: Path = ""
    output_datafile: Path = ""
    feff_file: list = field(factory=list)
    log_file: str = "test.csv"

    firstPass: bool = False
    multi_data_toggle: bool = False
    multi_data: list = field(factory=list)
    nComp: int = 1

    front: list = field(factory=list)

    data_path: Path = None
    output_path: Path = None
    log_path: Path = None

    pathOptimize: bool = False
    end: str = ".dat"

    def initialize_filepath(self, cycles: float = 0) -> None:
        """
        Initialize File Path
        """

        if self.nComp > 1:
            for i in range(self.nComp):
                self.front.append(self.base / self.feff_file[i])
        else:
            self.front = self.base / self.feff_file

        if self.multi_data_toggle:
            self.data_path = self.base / self.multi_data[cycles]
            # self.output_path = os.path.splitext(os.path.join(self.base, output_file))[
            #                        0] + "_" + str(i) + ".csv"
            self.output_path = Path(
                str((self.base / self.output_file).with_suffix("")) + f"_{cycles}.csv"
            )
            # self.log_path = os.path.splitext(
            #     copy.deepcopy(self.output_path))[0] + ".log"
            self.log_path = Path(str(self.output_path.with_suffix("")) + ".log")
        else:
            # self.data_path = os.path.join(self.base, csv_file)
            self.data_path = self.base / self.data_file
            self.output_path = self.base / self.output_file
            self.log_path = Path(str(self.output_path.with_suffix("")) + ".log")
            # self.log_path = os.path.splitext(
            #     copy.deepcopy(self.output_path))[0] + ".log"
        if self.pathOptimize:
            self.output_path = Path(
                str((self.base / self.output_file).with_suffix("")) + "_optimizes.csv"
            )
        #     self.output_path = os.path.splitext(os.path.join(self.base, output_file))[
        #                            0] + "_optimized.csv"

        self.initialize_outputs()

    def read_inputs(self, input_dicts: Dict) -> None:
        self.nComp = checkKey("num_compounds", input_dicts, 1)
        if self.nComp > 1:
            try:
                self.feff_file = list(input_dicts["feff_file"].split(","))
            except FileNotFoundError:
                print("Feff folder is not correct")
        else:
            self.feff_file = input_dicts["feff_file"]

        self.data_file = checkKey("data_file", input_dicts, "")
        self.output_file = checkKey("output_file", input_dicts, "exafs_neo_out.csv")
        self.log_file = checkKey("log_file", input_dicts, "exafs_neo.log")
        self.pathOptimize = checkKey("pathOptimize", input_dicts, False)

    def initialize_outputs(self) -> None:
        base_file = self.output_path.stem
        self.output_datafile = self.output_path.with_name(f"{base_file}_data.csv")

    def write_outputs(self, neoRunPars, bestFitPars) -> None:
        with open(self.output_path, "a") as f1:
            data_line = (
                f"{neoRunPars.currGen},{neoRunPars.tt},{bestFitPars.globBestVal}\n"
            )
            f1.writelines(data_line)

    def write_data_outputs(self, bestFitPars) -> None:
        globBestFit = bestFitPars.globBestInd.get()
        with open(self.output_datafile, "a") as f2:
            bestFit = globBestFit
            for path in bestFit:
                # f2.writerow(f"{i}")
                line = f"{path[0]},{path[1]},{path[2]},{path[3]}\n"
                f2.writelines(line)
            f2.write("#################################\n")


if __name__ == "__main__":
    # exafs_pars = EXAFSPars()
    inputs_pars = {
        "data_file": "path_files/Cu/cu_10k.xmu",
        "output_file": "tests/output.csv",
        "feff_file": "test/feff",
        "kmin": 0.95,
        "kmax": 9.775,
        "kweight": 3.0,
        "deltak": 0.05,
        "rbkg": 1.1,
        "bkgkw": 1.0,
        "bkgkmax": 15.0,
        "pathOptimize": False,
    }

    exafs_NeoPars = NeoFilePars()

    exafs_NeoPars.read_inputs(inputs_pars)
    exafs_NeoPars.initialize_filepath()
    print(exafs_NeoPars)
