from typing import Dict

import larch

try:
    from larch import Interpreter
    from larch.io import read_ascii
    from larch.xafs import autobk, feffdat, xftf
except:
    from larch_plugins.io import read_ascii
    from larch_plugins.xafs import autobk
    from larch_plugins.xafs import feffdat
    from larch_plugins.xafs import xftf

from attrs import define, field
import numpy as np
from scipy.integrate import simpson

from PhysicsModules.EXAFS.exafs_neo.pathrange import Pathrange_limits, read_pathrange_file

from PhysicsModules.EXAFS.exafs_neo.neoFilePars import NeoFilePars
from PhysicsModules.EXAFS.exafs_neo.utils import checkKey, time_call, Utils


@define
class NeoBestFit:
    currBestInd: float = None
    currBestVal: float = np.inf
    globBestInd: float = None
    globBestVal: float = np.inf

    bestDiff: float = np.inf
    bestE0: float = 0

    bestChir_magTotal: np.array = field(default=np.zeros(326))
    bestYTotal: np.array = field(default=np.zeros(326))


@define(kw_only=True, slots=True)
class NeoFixedPars:
    # Neo Fixed Parameters

    nPops: float = 100
    nGen: float = 100
    steadyState: bool = False

    selOpt: int = 0
    croOpt: int = 0

    DE: bool = False

    pathOptimize: bool = False
    pathOptimizePercent: float = 0.01
    pathOptimizeOnly: bool = False
    filter_percentage: float = 0.2

    solver: str = "GA"

    printGraph: bool = True
    debug_mode: bool = False

    def read_inputs(self, input_dicts):

        self.nPops = checkKey("nPops", input_dicts, 100)
        self.nGen = checkKey("nGen", input_dicts, 100)
        self.steadyState = checkKey("steadyState", input_dicts, False)

        self.selOpt = checkKey("selOpt", input_dicts, 0)
        self.croOpt = checkKey("croOpt", input_dicts, 0)

        self.DE = checkKey("DE", input_dicts, False)

        self.pathOptimize = checkKey("pathOptimize", input_dicts, False)
        self.pathOptimizePercent = checkKey("pathOptimize_percent", input_dicts, 0.01)
        self.pathOptimizeOnly = checkKey("pathOptimizeOnly", input_dicts, False)
        self.filter_percentage = checkKey("filter_percentage", input_dicts, 0.2)

        self.solver = checkKey("solver", input_dicts, "GA")

        self.printGraph = checkKey("printGraph", input_dicts, True)
        self.debug_mode = checkKey("debug_mode", input_dicts, False)

        if self.printGraph:
            # TODO: Implement this
            # self.fig = plt.figure()
            # self.ax = self.fig.add_subplot(111)
            pass


@define
class NeoRunPars:
    currGen: int = 1
    time: bool = False
    tt: float = 0
    secondHalf: bool = False
    diffCounter: int = 0
    cycles: int = 0
    currGen_st: float = 0
    currGen_tt: float = 0
    nGen: int = 100

    def start_gen(self):
        self.currGen_st = time_call()

    def end_gen(self, neo_population):
        self.currGen += 1
        self.currGen_tt = time_call() - self.currGen_st
        self.tt += self.currGen_tt
        if self.currGen == self.nGen // 2:
            neo_population.optimize_e0()
            self.secondHalf = True

        if self.currGen == self.nGen:
            neo_population.optimize_e0()

    def read_inputs(self, input_dicts):
        self.nGen = checkKey("nGen", input_dicts, 100)


@define
class NeoMutPars:
    """
    Neo Mutation Parameters
    """

    mutOpt: int = 1
    mutChance: float = field(default=0.3)
    mutChanceE0: float = field(default=0.3)
    nmut: int = 0

    def read_inputs(self, input_dicts):
        self.mutOpt = checkKey("mut_options", input_dicts, 1)
        self.mutChance = checkKey("mutChance", input_dicts, 0.3)
        self.mutChanceE0 = checkKey("mutChanceE0", input_dicts, 0.3)

    @mutChance.validator
    def check_mutchance(self, attribute, value):
        if value > 1.0 or value < 0.0:
            raise ValueError("mutChance should be between 0 and 1")


@define
class NeoCrossPars:
    croOpt: int = 0

    def read_inputs(self, input_dicts):
        self.croOpt = checkKey("croOpt", input_dicts, 0)


@define
class NeoSelPars:
    """
    Neo Selection Parameters
    """

    selOpt: int = 0
    nBestSample: float = 0.3
    nLuckSample: float = 0.2

    parents: list = field(factory=list)
    __nPop: int = 0
    nBest: int = 0
    nLuck: int = 0
    nCross: int = 0

    def read_inputs(self, input_dicts):
        self.selOpt = checkKey("selOpt", input_dicts, 0)
        self.nBestSample = checkKey("nBestSample", input_dicts, 0.3)
        self.nLuckSample = checkKey("nLuckSample", input_dicts, 0.2)
        self.__nPop = checkKey("nPops", input_dicts, 100)
        self.nBest = int(self.__nPop * self.nBestSample)
        self.nLuck = int(self.__nPop * self.nLuckSample)
        self.nCross = self.__nPop - self.nBest - self.nLuck


@define
class NeoSol:
    solOpt: int = 0

    def read_inputs(self, input_dicts):
        self.solOpt = checkKey("solOpt", input_dicts, 0)


class NeoPars:
    def __init__(self, verbose_lvl=5):
        """
        Wrapped all paras together
        """

        self.verbose_lvl = verbose_lvl
        self.fixedPars = NeoFixedPars()
        self.runPars = NeoRunPars()
        self.mutPars = NeoMutPars()
        self.crossPars = NeoCrossPars()
        self.selPars = NeoSelPars()
        self.exafsPars = EXAFSStaticPars()
        self.bestFitPars = NeoBestFit()
        self.neoFilePars = NeoFilePars()
        self.exafsRangePars = EXAFSPathRange()
        self.exafsPathPars = EXAFSPath()
        self.solPars = NeoSol()

    def read_inputs(self, input_dicts):
        self.fixedPars.read_inputs(input_dicts)
        self.runPars.read_inputs(input_dicts)
        self.exafsPars.read_inputs(input_dicts)
        self.neoFilePars.read_inputs(input_dicts)
        self.mutPars.read_inputs(input_dicts)
        self.selPars.read_inputs(input_dicts)
        self.solPars.read_inputs(input_dicts)
        self.crossPars.read_inputs(input_dicts)
        self.neoFilePars.initialize_filepath(cycles=0)

        self.exafsRangePars.read_inputs(input_dicts, self.exafsPars)
        self.exafsPathPars.read_inputs(self.neoFilePars, self.exafsPars)

        self.exafsPathPars.initialize()

    def output(self):
        self.neoFilePars.write_outputs(self.runPars, self.bestFitPars)
        self.neoFilePars.write_data_outputs(self.bestFitPars)

    def end_gen(self, neo_population):
        self.output()
        self.runPars.end_gen(neo_population)
        # self.exafsPathPars.end_gen()


@define
class EXAFSPathRange:
    pathrange_file: str = ""
    pathrange_pars: list = field(default=[])
    npath: int = 0

    # e0, for anything
    # rangeE0: np.array = field(default_factory=(np.linspace(-100, 100, 201) * 0.01))
    rangeE0: np.array = field(default=np.linspace(-100, 100, 201) * 0.01)
    # Large range B
    rangeE0_large: np.array = field(default=np.linspace(-600, 600, 1201) * 0.01)

    def read_inputs(self, input_dicts, exafs_pars):
        self.npath = exafs_pars.npath
        self.pathrange_file = checkKey("pathrange_file", input_dicts, None)
        # self.pathrange_pars = checkKey('pathrange_pars', input_dicts, [])
        if self.pathrange_file is None:
            for i in range(self.npath):
                self.pathrange_pars.append(Pathrange_limits(i))

        else:
            pathrange_file = read_pathrange_file(self.pathrange_file, self.npath)
            for i in range(self.npath):
                self.pathrange_pars.append(Pathrange_limits(i, pathrange_file[i, :]))

    def modify_pathrange(self, best_individual):
        for i in range(self.npath):
            path_bestfit = best_individual.get_path(i)
            self.pathrange_pars[i].mod_s02(path_bestfit[0])
            self.pathrange_pars[i].mod_sigma2(path_bestfit[2])
            self.pathrange_pars[i].mod_deltaR(path_bestfit[3])


@define(slots=True)
class EXAFSStaticPars:
    kmin: float = 0.95
    kmax: float = 9.775
    dk: float = 0.05
    kweight: float = 2.0

    rbkg: float = 0.0
    bkgkw: float = 1.0
    bkgkmax: float = 15.0

    small: float = 0.0
    big: float = 0.0
    mid: float = 0.0
    intervalK: list = field(factory=list)

    pathOptimize: bool = False
    pathOptimizePercent: float = 0.01
    pathOptimizeOnly: bool = False

    individual_paths: bool = False

    ncomp: int = 1
    pathrange: list = field(factory=list)
    npath: int = 0

    def calculate_pars(self):
        self.small = int(self.kmin / self.dk)
        self.big = int(self.kmax / self.dk)
        self.mid = int(self.big - self.small + 1)
        self.intervalK = np.linspace(self.small, self.big, self.mid)

    def read_inputs(self, input_dicts):
        self.kmin = checkKey("kmin", input_dicts, 0.95)
        self.kmax = checkKey("kmax", input_dicts, 9.775)
        self.dk = checkKey("deltak", input_dicts, 0.05)
        self.kweight = checkKey("kweight", input_dicts, 2.0)

        self.rbkg = checkKey("rbkg", input_dicts, 0.0)
        self.bkgkw = checkKey("bkgkw", input_dicts, 1.0)
        self.bkgkmax = checkKey("bkgkmax", input_dicts, 15.0)

        self.pathrange = checkKey("pathrange", input_dicts, None)
        self.ncomp = checkKey("num_compounds", input_dicts, 1)

        if self.ncomp > 1:
            self.individual_paths = True
            self.npath = sum(len(comp) for comp in self.pathrange)
        else:
            self.individual_paths = self.individual_paths
            self.npath = len(self.pathrange)

        self.calculate_pars()

        self.pathOptimize = checkKey("pathOptimize", input_dicts, False)
        self.pathOptimizePercent = checkKey("pathOptimize_percent", input_dicts, 0.01)
        self.pathOptimizeOnly = checkKey("pathOptimizeOnly", input_dicts, False)


@define
class EXAFSPath:
    mylarch: str = Interpreter()
    g: larch.symboltable.Group = None
    best: larch.symboltable.Group = None
    sumgroup: larch.symboltable.Group = None
    exp: list = field(factory=list)
    pathname: list = field(factory=list)
    ncomp: int = 1
    individual_paths: bool = False
    path_lists: list = field(factory=list)
    pathDictionary: dict = field(factory=dict)
    end: str = None
    front: list = field(factory=list)
    npath: int = None
    exafs_static_pars: EXAFSStaticPars = None
    exafs_file_pars: NeoFilePars = None

    def read_inputs(self, exafs_filepars, exafs_static_pars: EXAFSStaticPars):
        self.exafs_file_pars = exafs_filepars
        self.exafs_static_pars = exafs_static_pars
        self.g = read_ascii(str(exafs_filepars.data_path))
        self.best = read_ascii(str(exafs_filepars.data_path))
        self.sumgroup = read_ascii(str(exafs_filepars.data_path))
        self.ncomp = exafs_filepars.nComp
        self.individual_paths = exafs_static_pars.individual_paths
        self.npath = exafs_static_pars.npath
        self.front = exafs_filepars.front
        self.end = exafs_filepars.end
        self.path_lists = exafs_static_pars.pathrange

    def __initialize_group(self):
        try:
            self.g.k
            self.g.chi
        except AttributeError:
            autobk(
                self.g,
                rbkg=self.exafs_static_pars.rbkg,
                kweight=self.exafs_static_pars.bkgkw,
                kmax=self.exafs_static_pars.bkgkmax,
                _larch=self.mylarch,
            )
            autobk(self.best, rbkg=self.exafs_static_pars.rbkg, _larch=self.mylarch)
            autobk(self.sumgroup, rbkg=self.exafs_static_pars.rbkg, _larch=self.mylarch)

    def __initialize_paths(self):
        """
        Load paths:
            Initialize paths in various files.
        """
        self.pathname = []
        if self.ncomp > 1:
            if self.individual_paths:
                for i in range(self.ncomp):
                    comp_path = len(self.path_lists[i])
                    for j in range(comp_path):
                        filename = (
                            str(self.front[i])
                            + str(self.path_lists[i][j]).zfill(4)
                            + self.end
                        )

                        pathName = "Comp" + str(i) + "Path" + self.path_lists[i][j]
                        self.pathname.append(pathName)
                        self.pathDictionary.update(
                            {pathName: feffdat.feffpath(filename, _larch=self.mylarch)}
                        )
        else:
            if self.individual_paths:
                for i in range(self.npath):
                    filename = (
                        str(self.front) + str(self.path_lists[i]).zfill(4) + self.end
                    )
                    pathName = "Path" + self.path_lists[i]
                    self.pathname.append(pathName)
                    self.pathDictionary.update(
                        {pathName: feffdat.feffpath(filename, _larch=self.mylarch)}
                    )

            else:
                for i in range(1, self.npath + 1):
                    filename = str(self.front) + str(i).zfill(4) + self.end
                    pathName = f"Path{i}"
                    self.pathname.append(pathName)
                    self.pathDictionary.update(
                        {pathName: feffdat.feffpath(filename, _larch=self.mylarch)}
                    )

    def __initialize_ftf(self):
        xftf(
            self.g.k,
            self.g.chi,
            kmin=self.exafs_static_pars.kmin,
            kmax=self.exafs_static_pars.kmax,
            dk=4,
            window="hanning",
            kweight=self.exafs_static_pars.kweight,
            group=self.g,
            _larch=self.mylarch,
        )
        xftf(
            self.best.k,
            self.best.chi,
            kmin=self.exafs_static_pars.kmin,
            kmax=self.exafs_static_pars.kmax,
            dk=4,
            window="hanning",
            kweight=self.exafs_static_pars.kweight,
            group=self.best,
            _larch=self.mylarch,
        )
        xftf(
            self.sumgroup.k,
            self.sumgroup.chi,
            kmin=self.exafs_static_pars.kmin,
            kmax=self.exafs_static_pars.kmax,
            dk=4,
            window="hanning",
            kweight=self.exafs_static_pars.kweight,
            group=self.sumgroup,
            _larch=self.mylarch,
        )

        self.exp = self.g.chi

    @staticmethod
    def path_optimization_multicomp(path_lists, new_path_ind, ncomp):
        """Path optimizations for multiple components

        Args:
            path_lists (lists): Path lists
            new_path_ind (lists): Optimize path index
            ncomp (int): Number of components

        Returns:
            (list): Optimized path lists
        """
        flat_path = [element for sublist in path_lists for element in sublist]

        cumsum_path = np.cumsum(Utils.get_sum_path(path_lists, ncomp))

        new_path_list = []
        temp_path = []

        k = 0
        for i, path in enumerate(flat_path):
            if i in new_path_ind:
                temp_path.append(path)
            if i + 1 == cumsum_path[k]:
                new_path_list.append(temp_path)
                k += 1
                temp_path = []

        return new_path_list

    def path_optimize(
        self, best_individual, exafs_static_pars: EXAFSStaticPars, path_optimize_percent
    ) -> Dict:
        arr = np.asarray(best_individual.get())
        total = 0
        # Contribution of each path
        contrib = []
        for i, (key, value) in enumerate(self.pathDictionary.items()):
            path = feffdat.feffpath(
                value.filename,
                s02=str(arr[i, 0]),
                e0=str(arr[i, 1]),
                sigma2=str(arr[i, 2]),
                deltar=str(arr[i, 3]),
                _larch=self.mylarch,
            )
            feffdat.path2chi(path, _larch=self.mylarch)

            xftf(
                path.k,
                path.chi,
                kmin=exafs_static_pars.kmin,
                kmax=exafs_static_pars.kmax,
                dk=4,
                window="hanning",
                kweight=exafs_static_pars.kweight,
                group=path,
                _larch=self.mylarch,
            )

            area = simpson(path.chir_mag, path.r)
            total += area
            contrib.append(area)
        # Contrib percentage
        contrib_p = [i / total for i in contrib]
        new_path_ind = (
            np.argwhere(np.array(contrib_p) >= path_optimize_percent)
        ).flatten()

        new_path_dict = {
            "path_optimize_percent": path_optimize_percent,
            "original_paths": self.path_lists,
        }
        optimized_path = {}
        if self.ncomp > 1:
            new_path = EXAFSPath.path_optimization_multicomp(
                self.path_lists, new_path_ind, self.ncomp
            )

            for i, comp_list in enumerate(new_path):
                this_comp_str = f"Comp{i}"
                for j, path in enumerate(comp_list):
                    path_ind = i * len(comp_list) + j
                    optimized_path[f"{this_comp_str}Path{path}"] = float(
                        contrib_p[path_ind]
                    )
        else:
            new_path = [self.path_lists[i] for i in new_path_ind if i in new_path_ind]
            for i, path in enumerate(self.path_lists):
                if i in new_path_ind:
                    optimized_path[str(path)] = float(contrib_p[i])
        new_path_dict["optimize_paths"] = new_path
        new_path_dict["optimize_path_lists"] = optimized_path
        return new_path_dict

    def initialize(self):
        self.__initialize_group()
        self.__initialize_paths()
        self.__initialize_ftf()


if __name__ == "__main__":
    # exafs_pars = EXAFSPars()
    inputs_pars = {
        "data_file": "../path_files/Cu/cu_10k.xmu",
        "output_file": "",
        "feff_file": "../path_files/Cu/path_75/feff",
        "kmin": 0.95,
        "kmax": 9.775,
        "kweight": 3.0,
        "pathrange": [1, 2, 3, 4, 5],
        "deltak": 0.05,
        "rbkg": 1.1,
        "bkgkw": 1.0,
        "bkgkmax": 15.0,
    }

    # exafs_NeoPars = NeoFilePars()
    # exafs_NeoPars.read_inputs(inputs_pars)
    # exafs_NeoPars.initialize_filepath()
    # print(exafs_NeoPars)
    exafs_NeoPars = NeoPars()
    exafs_NeoPars.read_inputs(inputs_pars)

    print(exafs_NeoPars.exafsRangePars)
    # print(exafs_NeoPars.exafsPathPars)
