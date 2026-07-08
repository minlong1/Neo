"""
Authors    Miu Lun(Andy) Lau*, Jeffrey Terry, Min Long
Email      andylau@u.boisestate.edu, jterry@agni.phys.iit.edu, minlong@boisestate.edu
Version    0.2.0
Date       July 4, 2021

EXAFS Analysis Wrapper Functions, contains Latex, Igor Pro
"""

import fnmatch
import os
import re
from pathlib import Path
from typing import List, Optional, Dict, Union, TextIO
import matplotlib.pyplot as plt
import numpy as np
from larch.xafs import xftf
from matplotlib.figure import Figure
from scipy.integrate import simpson
from PhysicsModules.EXAFS.exafs_neo.analysis import larch_score
from PhysicsModules.EXAFS.exafs_neo.gui.Misc_Function import natural_keys


def calculate_occurrences(folder_name: str, paths: List[str], limits=20.0) -> List[int]:
    """
    Calculate Occurrences plot
    @param str folder_name: Folder name
    @param list paths: path list
    @param float limits: limit parameters. Defaults to 20.
    @return:
    """
    #  calculate the number of log for each path optimization
    files = []
    for r, d, f in os.walk(folder_name):
        f.sort(key=natural_keys)
        for file in f:
            if fnmatch.fnmatch(file, "*.log"):
                files.append(os.path.join(r, file))

    occ_list = np.zeros(len(paths))
    for i in range(len(files)):
        with open(files[i]) as f:
            for line in f:
                if "New Paths:" in line:
                    temp_str = []
                    list_str = line[12:-2].split(",")
                    for j in range(len(list_str)):
                        temp_str.append(int(list_str[j].replace("'", "")))
                    for _, cur_path in enumerate(temp_str):
                        for j, og_path in enumerate(folder_name):
                            if cur_path == og_path:
                                occ_list[j] += 1

    return list(occ_list)


def plot_occ_list(folder_name: str, limits: int, paths: List, fig_gui=None) -> None:
    """Plot Occurrences list

    Args:
        folder_name (str): folder name
        limits (float): limit parameters, default is 20
        paths (list): _description_
        fig_gui (fig, optional): _description_. Defaults to None.
    """
    occ_list = calculate_occurrences(folder_name, paths, limits)
    if fig_gui is None:
        plt.figure(figsize=(8, 5))
        plt.xticks(paths)
        plt.bar(paths, occ_list)
    else:
        ax = fig_gui.add_subplot(111)
        ax.set_xticks(paths)
        # ax.set_xticklabels(self.label,rotation=70)
        ax.bar(paths, occ_list)
        fig_gui.tight_layout()


class EXAFS_Analysis:
    """
    Set up analysis of EXAFS dataset

    Inputs:
        paths (list): list containing paths numbers
        dirs (str): str pointing the files locations
        params (dics): dicts contain all parameters
    """

    def __init__(self, paths: List[List[int]], dirs: str, params: Dict) -> None:
        self.ind_export_paths = None
        self.paths = paths
        # Generate Label for plotting
        self.flat_paths = larch_score.flatten_2d_list(self.paths)
        self.label = larch_score.generate_labels(self.flat_paths)[0]
        self.optimize_path = None
        self.params = params
        self.dirs = dirs
        self.num_paths = len(self.flat_paths)
        self.num_params = 3 * self.num_paths + 1
        self.return_str = ""
        self.params["optimize"] = params.get("optimize", False)
        if "series" not in params:
            self.series = False
        else:
            self.series = self.params["series"]
            self.series_index = self.params["series_index"]

    def read_optimize_paths(self, files_opt) -> None:
        """
        Read optimize paths

        Args:
            files_opt (file): optimized paths
        """
        optimize_path = []
        for i in range(len(files_opt)):
            with open(files_opt[i]) as fp:
                for j, line in enumerate(fp):
                    if j == 1:
                        # clean up the newline symbol
                        clean = line.replace("]\n", "")
                        clean = clean.replace("[", "")
                        clean_arr = np.fromstring(clean, dtype=int, sep="  ")

                        optimize_path.append(clean_arr)

        self.optimize_path = optimize_path

    def read_result_files(
        self, series: bool = False, series_index: int = None, verbose: bool = False
    ):
        num_path = self.num_paths
        full_mat = []
        best_full_mat = None
        files = []
        folder = self.dirs

        if self.params["optimize"]:
            files_opt = []
            files_opt_data = []

        if os.path.isdir(folder):
            search_string = "*_*_" + str(series_index) + "_data.csv"

            f = os.listdir(folder)
            f.sort(key=natural_keys)
            for file in f:
                # various checks in place
                if not series:
                    if fnmatch.fnmatch(file, "*_data.csv"):
                        files.append(os.path.join(folder, file))
                elif series:
                    if fnmatch.fnmatch(file, search_string):
                        files.append(os.path.join(folder, file))
                elif self.params["optimize"]:
                    if re.search(r"test_\d+_optimized.csv", file):
                        files_opt.append(os.path.join(folder, file))
                    if re.search(r"test_\d+_optimized_data.csv", file):
                        files_opt_data.append(os.path.join(folder, file))
        else:
            self.params["individual"] = True
            files.append(folder)
        if verbose:
            print(folder)
        files.sort(key=natural_keys)
        if self.params["optimize"]:
            files_opt.sort(key=natural_keys)
            files_opt_data.sort(key=natural_keys)
            self.read_optimize_paths(files_opt)
        # Loop through all the files.
        for i in range(len(files)):
            file = files[i]
            try:
                # Check if the file path exists
                os.path.exists(file)
                # Read the resulting csv
                gen_csv = np.genfromtxt(file, delimiter=",")
                gen_csv_unflatten = gen_csv.reshape((-1, 4 * num_path))

                gen_csv_best_unflatten = gen_csv[-num_path::].reshape(
                    (-1, 4 * num_path)
                )
                try:
                    self.params["individual"]
                except KeyError:
                    self.params["individual"] = False

                if self.params["individual"]:
                    full_mat = gen_csv_unflatten
                    best_full_mat = gen_csv_best_unflatten
                else:
                    if i == 0:
                        full_mat = gen_csv_unflatten
                        best_full_mat = gen_csv_best_unflatten
                    best_full_mat = np.vstack((best_full_mat, gen_csv_best_unflatten))
                    full_mat = np.vstack((full_mat, gen_csv_unflatten))
            except OSError:
                print(" " + str(i) + " Missing")
                pass

        if best_full_mat is None:
            raise FileNotFoundError(
                f"No *_data.csv result files found under {folder!r} "
                f"(dir exists: {os.path.isdir(folder)}, matched {len(files)} file(s))"
            )
        return full_mat, best_full_mat

    def extract_data(
        self,
        data: List[List[float]] = None,
        plot: bool = False,
        average: bool = False,
        verbose: bool = True,
    ):
        """
        Extract the data from a list of  dataset and construct a distribution
        """

        if data is not None:
            full_mat = data["full_mat"]
            bestfit_full_mat = data["bestfit_full_mat"]
        else:
            full_mat, bestfit_full_mat = self.read_result_files(
                self.series, self.series_index, verbose=verbose
            )
        self.full_mat = full_mat
        self.bestfit_full_mat = bestfit_full_mat
        self.create_label()
        self.construct_err_mat(full_mat, plot)
        self.construct_bestfit_mat(self.bestfit_full_mat, average=average)

    def larch_init(self) -> None:
        """
        Initialize larch helper function in the background
        """
        csv_path = os.path.join(self.params["base"], self.params["CSV"])
        exp, g, params, mylarch = larch_score.larch_init(csv_path, self.params)

        self.g = g
        self.exp = exp
        self.small = params.get("SMALL")
        self.big = params.get("BIG")
        self.intervalK = params.get("intervalK")
        self.kweight = params.get("kweight")
        self.mylarch = mylarch

    def larch_score(self, verbose: bool = True):
        """
        Calculate fitness score based on Chi and ChiR

        Args:
            verbose (bool, optional): if verbose output. Defaults to True.
        """

        (
            self.path,
            self.yTotal,
            self.best,
            self.loss,
            self.best_Fit_r,
            self.og_arr_r,
            arr,
        ) = larch_score.fitness(
            self.exp,
            self.bestFit_mat,
            self.paths,
            self.params,
            self.mylarch,
            self.g,
            verbose=verbose,
        )
        self.return_str += arr

        self.chir2 = (self.loss / (len(self.intervalK) - self.num_params)).round(6)

        if verbose:
            print(f"Fitness Score (Chi2): {np.round(self.loss, 6)}")
            print(f"Fitness Score (ChiR2): {np.round(self.chir2, 6)}")

    def plot(
        self,
        title: Optional[str] = "Temp",
        fig_gui: Figure = None,
        show: bool = False,
        save_file: Union[str, Path] = None,
    ) -> None:

        SMALL = self.small
        BIG = self.big

        self.best.k = self.path.k
        self.best.chi = self.yTotal
        xftf(
            self.best.k,
            self.best.chi,
            kmin=self.params["Kmin"],
            kmax=self.params["Kmax"],
            dk=4,
            window="hanning",
            kweight=self.kweight,
            group=self.best,
        )

        if fig_gui is None:
            if save_file is not None:
                show = True
            if show:
                fig, ax = plt.subplots(1, 2, figsize=(15, 5))

                ax[0].plot(
                    self.g.k,
                    self.g.chi * self.g.k**self.kweight,
                    "b--",
                    label="Experiment Data",
                )
                ax[0].plot(
                    self.path.k[SMALL:BIG],
                    self.yTotal[SMALL:BIG] * self.path.k[SMALL:BIG] ** self.kweight,
                    "r-",
                    label="Genetic Algorithm",
                )
                ax[0].set_xlim((0, self.params["Kmax"]))
                ax[0].legend()
                ax[0].set_title(title + " K Space")
                ax[0].set_xlabel(r"$k$ (Å$^{-1}$)")
                ax[0].set_ylabel(
                    r"k$^{"
                    + str(int(self.kweight))
                    + r"}\chi$(k) (Å$^{"
                    + str(int(-1 * self.kweight))
                    + "}$)"
                )

                ax[1].plot(self.g.r, self.g.chir_mag, "b--", label="Experiment Data")
                ax[1].plot(
                    self.best.r, self.best.chir_mag, "r-", label="Genetic Algorithm"
                )
                ax[1].set_title(title + " R Space")
                ax[1].legend()
                ax[1].set_xlabel(r"$R$ (Å)")
                negative_weight = str(-1 * (int(self.kweight) + 1))
                ax[1].set_ylabel(r"FT[$\chi$(r)](Å$^{" + negative_weight + R"}$)")
                if save_file is not None:
                    fig.savefig(save_file, dpi=500)
        else:
            ax = fig_gui.add_subplot(121)

            ax.plot(
                self.g.k,
                self.g.chi * self.g.k**self.kweight,
                "b--",
                label="Experiment Data",
            )
            ax.plot(
                self.path.k[SMALL:BIG],
                self.yTotal[SMALL:BIG] * self.path.k[SMALL:BIG] ** self.kweight,
                "r-",
                label="Genetic Algorithm",
            )
            ax.legend()
            ax.set_title(title + " K Space")

            ax = fig_gui.add_subplot(122)

            ax.plot(self.g.r, self.g.chir_mag, "b--", label="Experiment Data")
            ax.plot(self.best.r, self.best.chir_mag, "r-", label="Genetic Algorithm")
            ax.set_title(title + " R Space")
            ax.legend()
            fig_gui.tight_layout()

    def construct_latex_table(self, print_table: bool = False) -> None:
        """
        Construct latex table

        Args:
            print_table (bool, optional): print the table. Defaults to False.


        Todo:
            Change this from printout to files instead.
        """

        err_full = larch_score.construct_full_err(self.err)
        nleg_arr, label_arr, latex_table_str = larch_score.latex_table(
            self.paths, self.best_Fit_r, err_full, self.chir2, self.og_arr_r
        )

        self.nleg_arr = nleg_arr
        self.label_arr = label_arr
        self.latex_table_str = latex_table_str
        if print_table:
            print(self.latex_table_str)

    def individual_fit(self, plot: bool = False, fig_gui: Figure = None) -> None:
        """
        Perform fitness calculation separately for each path for r space

        Args:
            plot (bool, optional): plotting. Defaults to False.
            fig_gui (_type_, optional): supplement gui for plotting in other windows. Defaults to None.
        """

        # path,yTotal,best,loss,export_paths=larch_score.fitness_individual(self.exp,self.best_Fit,self.paths,self.params,export=True,plot=True)
        path, yTotal, best, loss, export_paths = larch_score.fitness_individual(
            self.exp,
            self.bestFit_mat,
            self.paths,
            self.params,
            mylarch=self.mylarch,
            g = self.g,
            export=True,
            plot=plot,
            fig_gui=fig_gui,
        )
        self.ind_export_paths = export_paths

    def write_latex_csv(self, file_name: str = "") -> None:
        """
        Write latex table to csv file
        """

        # check if attributes exists first.
        if hasattr("self", "latex_table_str"):
            print("Attribute Exists")
        else:
            self.construct_latex_table()
        if file_name is not None:
            with open(file_name, mode="w", newline="", encoding="utf-8") as write_file:
                write_file.write(self.latex_table_str)

    def get_data(
        self,
    ) -> tuple[tuple[List[float], List[float], List[float], List[float]], List[float]]:
        """
        Get the data for plotting

        Returns:
            tuple: (data_x, data_y, model_x, model_y),error_full
        """
        small = self.small
        big = self.big
        data_x = self.g.k[small:big]
        data_y = self.g.chi[small:big] * self.g.k[small:big] ** self.kweight
        model_x = self.path.k[small:big]
        model_y = self.yTotal[small:big] * self.path.k[small:big] ** self.kweight

        error_full = larch_score.construct_full_err(self.err)
        return (data_x, data_y, model_x, model_y), error_full

    def export_files(
        self, header: str = "test", dirs: Union[str, Path] = "", igor_true: bool = False
    ) -> None:
        """Export files to their corresponding spots:

        Args:
            header (str, optional): header file. Defaults to 'test'.
            dirs (str, optional): path location for files. Defaults to ''.
            igor_true (bool, optional): if igor plotting. Defaults to False.

        Todo:
            Removed one of the files since it is plotted directly using latex
            function

        """

        file_name_k = os.path.join(dirs, "bestFit_" + header + ".csv")
        file_name_ind = os.path.join(dirs, "Individual_Fit_" + header + ".csv")
        file_name_latex = os.path.join(dirs, "Latex_" + header + ".txt")

        SMALL = self.small
        BIG = self.big
        self.header = header
        larch_score.write_bestFit_csv(
            self.g.k,
            self.g.chi,
            self.path.k,
            self.yTotal,
            SMALL,
            BIG,
            name=file_name_k,
            header_base=header,
        )
        # larch_score.write_result_n_err(full_mat_var,err,name=file_name_best,header_base=header)
        larch_score.write_individual_csv(
            self.g.k,
            self.g.chi,
            self.path.k,
            self.yTotal,
            SMALL,
            BIG,
            self.ind_export_paths,
            self.paths,
            name=file_name_ind,
            header_base=header,
        )

        self.write_latex_csv(file_name_latex)
        if igor_true:
            self.export_igor_individual()

    def construct_bestfit_mat(
        self, bestfit_full_mat: List[float], average: bool
    ) -> None:
        if average:
            bestfit_full_mat_var = np.delete(bestfit_full_mat, np.s_[5::4], axis=1)
            best_fit = np.mean(bestfit_full_mat_var, axis=0)
        else:
            chi2_arr = []
            for this_instance in bestfit_full_mat:
                this_instance_unflatten = this_instance.reshape(-1, 4)
                _, _, _, loss, _, _, _ = larch_score.fitness(
                    self.exp,
                    this_instance_unflatten,
                    self.paths,
                    self.params,
                    self.mylarch,
                    self.g,
                    verbose=False,
                )
                chi2_arr.append(loss)
                print(loss)
                print(this_instance_unflatten)
            # [print(f"\t File_{str(i).zfill(3)}: {val}") for i,val in enumerate(chi2_arr)]
            min_chi2 = np.argmin(chi2_arr)
            best_fit = bestfit_full_mat[min_chi2]
            best_fit_mat = best_fit.reshape(-1, 4)
        # TODO: Individuals current is wrong
        # if self.params["individual"]:
        #     best_fit = bestfit_full_mat.reshape(-1, 4).round(6)
        # else:
        #     best_fit_mat = self.full_mat[-1].reshape(-1, 4).round(6)

        self.bestFit = best_fit
        self.bestFit_mat = best_fit_mat

    def create_label(self) -> None:
        label = larch_score.generate_labels(self.paths)[0]
        self.label = label

    def construct_err_mat(
        self, full_mat: Union[List[float], np.ndarray], plot: bool = False
    ) -> None:
        """
        Construct the error_matrix .

        full_mat (np.array): flatten array of all iterations and samples

        """

        # Delete the duplicate columns for e0.
        full_mat_var = np.delete(full_mat, np.s_[5::4], axis=1)
        # bestfit_full_mat_var = np.delete(bestfit_full_mat, np.s_[5::4], axis=1)
        # Calculate the covariance matrix and extract the diagonal for error calculation
        full_mat_var_cov = np.cov(full_mat_var.T)
        full_mat_diag = np.diag(full_mat_var_cov)
        err = np.sqrt(full_mat_diag)
        self.err = err
        self.full_mat_diag = full_mat_diag

        if plot:
            plt.figure(figsize=(8, 5))
            plt.xticks(np.arange(len(full_mat_diag)), self.label, rotation=70)
            plt.bar(np.arange(len(full_mat_diag)), np.sqrt(full_mat_diag))

    def plot_error(self, fig_gui=None):
        if fig_gui is None:
            plt.figure(figsize=(8, 5))
            # plt.xticks(np.arange(len(self.full_mat_diag)),self.label,rotation=70);
            plt.bar(np.arange(len(self.full_mat_diag)), np.sqrt(self.full_mat_diag))
        else:
            ax = fig_gui.add_subplot(111)
            # ax.xticks
            ax.set_xticks(np.arange(len(self.full_mat_diag)))
            # ax.set_xticklabels(self.label,rotation=70)
            ax.bar(np.arange(len(self.full_mat_diag)), np.sqrt(self.full_mat_diag))
            fig_gui.tight_layout()

    def paths_optimizations(self, number=0.01, verbose=False):
        """
        Paths optimizations using simpsons area calculation.
        The calculation are used to perform

        Args:
            number (float, optional): Optimization cutoff percentage. Defaults to 0.01.
            verbose (bool, optional): verbosity. Defaults to False.
        """
        total = 0
        total_area = 0
        contrib = []
        contrib_area = []
        for i in range(self.num_paths):
            self.best.k = self.ind_export_paths[2 * i, :]
            self.best.chi = self.ind_export_paths[2 * i + 1, :]

            xftf(
                self.best.k,
                self.best.chi,
                kmin=self.params["Kmin"],
                kmax=self.params["Kmax"],
                dk=4,
                window="hanning",
                kweight=self.params["kweight"],
                group=self.best,
                _larch=self.mylarch,
            )

            total += np.linalg.norm(self.best.chir_mag)
            contrib.append(np.linalg.norm(self.best.chir_mag))
            contrib_area.append(simpson(self.best.chir_mag, x=self.best.r))
            total_area += simpson(self.best.chir_mag, x=self.best.r)
        contrib_p = [i / total for i in contrib]
        contrib_ap = [i / total_area for i in contrib_area]
        if verbose:
            print("Paths, Contrib Percentage (2-Norm), Contrib Percentage (Area)")
            for i in range(len(self.paths)):
                print(i + 1, contrib_p[i].round(3), contrib_ap[i].round(3))
        new_path = (np.argwhere(np.array(contrib_ap) >= number)).flatten() + 1
        print("New Paths:")
        print(new_path)
        plt.bar(np.arange(self.num_paths), height=contrib_ap)
        plt.xticks(np.arange(self.num_paths), self.flat_paths)

    def export_igor_individual(self, file_paths: str = "export_ind.ipf") -> None:
        """
        Export files in igor plotting for individuals, must be ran after
        the individual methods

        Args:
            file_paths (str, optional): location for igor plot script. Defaults to 'export_ind.ipf'.
        """
        # Displace all data
        f = open(file_paths, "w")
        f.write("•Display data_" + self.header + "_chi2 vs data_" + self.header + "_k;")
        f.write("\n")
        f.write(
            "•AppendToGraph fit_" + self.header + "_chi2 vs fit_" + self.header + "_k;"
        )
        f.write("\n")
        full_paths = larch_score.flatten_2d_list(self.paths)
        for i in range(len(full_paths)):
            f.write(
                "•AppendToGraph path_"
                + str(full_paths[i])
                + "_"
                + self.header
                + "_chi2 vs path_"
                + str(full_paths[i])
                + "_"
                + self.header
                + "_k;"
            )
            f.write("\n")

        f.write("•SetAxis bottom *,11")
        f.write("\n")

        ## Offset
        # offset first two is designated number
        # Todo:
        # Need to redesign this later!

        if len(full_paths) < 3:
            if len(full_paths) == 1:
                f.write(
                    "•ModifyGraph offset(path_"
                    + str(full_paths[0])
                    + "_"
                    + self.header
                    + "_chi2)={0,5}"
                )
                f.write("\n")
            elif len(full_paths) == 2:
                f.write(
                    "•ModifyGraph offset(path_"
                    + str(full_paths[0])
                    + "_"
                    + self.header
                    + "_chi2)={0,5}"
                )
                f.write("\n")
                f.write(
                    "•ModifyGraph offset(path_"
                    + str(full_paths[1])
                    + "_"
                    + self.header
                    + "_chi2)={0,10}"
                )
                f.write("\n")
        else:
            f.write(
                "•ModifyGraph offset(path_"
                + str(full_paths[0])
                + "_"
                + self.header
                + "_chi2)={0,5}"
            )
            f.write("\n")
            f.write(
                "•ModifyGraph offset(path_"
                + str(full_paths[1])
                + "_"
                + self.header
                + "_chi2)={0,10}"
            )
            f.write("\n")
            f.write(
                "•ModifyGraph offset(path_"
                + str(full_paths[2])
                + "_"
                + self.header
                + "_chi2)={0,12.5}"
            )
            f.write("\n")

        # offset the rest
        for i in range(len(full_paths) - 3):
            curr_paths = str(full_paths[i + 3])
            f.write(
                "•ModifyGraph offset(path_"
                + curr_paths
                + "_"
                + self.header
                + "_chi2)={0,"
                + str(15 + i)
                + "}"
            )
            f.write("\n")
        f.write(r'•Label left "k\\S2\\M \u03c7(k) (Å\\S-2\\M)";DelayUpdate')
        f.write("\n")
        f.write(r'•Label bottom "k (Å\\S-1\\M)"')
        f.write("\n")
        f.write("•ModifyGraph lsize(fit_" + self.header + "_chi2)=2")
        f.write("\n")
        for i in range(len(full_paths)):
            f.write(
                "•ModifyGraph lsize(path_"
                + str(full_paths[i])
                + "_"
                + self.header
                + "_chi2)=2"
            )
            f.write("\n")
        ## Legend
        # \r  - new line
        self.adjust_color(f)
        self.create_legend(f)
        f.write("•ModifyGraph mode(data_" + self.header + "_chi2)=3")
        f.write("\n")

    def create_legend(self, f: TextIO) -> None:
        """Create legend

        Args:
            f (writer): writer for igor file
        """
        legend_header = (
            r'•Legend/C/N=text0/J "Test\rk\\S2\\M\u03c7(k)\rTest_Detail\r\r\\s('
        )
        legend_1 = (
            "data_"
            + self.header
            + r"_chi2) Data\r\\s(fit_"
            + self.header
            + r'_chi2) Fit";'
        )
        legend = legend_header + legend_1
        f.write(legend)
        f.write("\n")
        full_paths = larch_score.flatten_2d_list(self.paths)
        for i in range(len(full_paths)):
            if int(self.nleg_arr[i]) > 2:
                addition = " MS "
            else:
                addition = " "
            paths_arr = (
                str(full_paths[i])
                + "_"
                + self.header
                + "_chi2) Path "
                + str(full_paths[i])
                + addition
                + self.label_arr[i]
                + r'";DelayUpdate'
            )
            f.write('•AppendText/N=text0 "\\s(path_' + paths_arr)
            f.write("\n")

    def adjust_color(self, f: TextIO, color_map=plt.cm.jet_r) -> None:
        """
        Adjust the color using jet reverse color bar
        """
        full_paths = larch_score.flatten_2d_list(self.paths)
        x = np.linspace(0, 1, len(full_paths))
        color = [color_map(i) for i in x]
        test = 65535
        for i in range(len(color)):
            color[i] = (
                int(test * color[i][0]),
                int(test * color[i][1]),
                int(test * color[i][2]),
            )
        # Change to X and Y for data
        f.write("•ModifyGraph rgb(fit_" + self.header + "_chi2)=(0,0,0)")
        f.write("\n")
        for i in range(len(color)):
            f.write(
                "•ModifyGraph rgb(path_"
                + str(full_paths[i])
                + "_"
                + self.header
                + "_chi2)="
                + str(color[i])
            )
            f.write("\n")

    def stacks_plot(self) -> None:
        """Generate stack plot"""
        self.best.k = self.ind_export_paths[0, :]
        self.best.chi = self.ind_export_paths[1, :]
        xftf(
            self.best.k,
            self.best.chi,
            kmin=self.params["Kmin"],
            kmax=self.params["Kmax"],
            dk=4,
            window="hanning",
            kweight=self.params["kweight"],
            group=self.best,
            _larch=self.mylarch,
        )
        y_arr = np.zeros((self.num_paths, len(self.best.r)))
        y_tot = np.zeros(len(self.best.r))
        # Get all the data
        for i in range(self.num_paths):
            self.best.k = self.ind_export_paths[2 * i, :]
            self.best.chi = self.ind_export_paths[2 * i + 1, :]

            xftf(
                self.best.k,
                self.best.chi,
                kmin=self.params["Kmin"],
                kmax=self.params["Kmax"],
                dk=4,
                window="hanning",
                kweight=self.params["kweight"],
                group=self.best,
                _larch=self.mylarch,
            )
            y_arr[i, :] = -self.best.chir_mag
            y_tot += self.best.chir_mag
        x = self.best.r

        plt.rc("font", size=11)
        rc = {"font.family": "serif", "mathtext.fontset": "stix"}
        plt.rcParams.update(rc)
        plt.rcParams["font.serif"] = ["Times New Roman"] + plt.rcParams["font.serif"]

        figsize = (10, 7)
        linewidth = 1.0
        edgecolor = "black"
        capsize = 5
        spacing = 0.2
        color = "royalblue"
        label_fontsize = 14

        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=figsize)
        ax.set_xlabel(r"R, uncorrected ($\AA$)")
        ax.set_ylabel(r"|$\chi$(R)| ($\AA^{3}$)")

        ax.plot(x, y_tot, "k", linestyle="solid", linewidth=linewidth, label="Fits")
        ax.stackplot(x, y_arr, labels=np.arange(1, self.num_paths + 1))

        # ax.set_xlim([1,9])
        ax.yaxis.grid(True, linestyle="--", which="major", color="grey", alpha=0.25)
        ax.xaxis.grid(True, linestyle="--", which="major", color="grey", alpha=0.25)
        ax.tick_params(which="both", direction="in")

        ax.legend(loc="lower right", fontsize=label_fontsize)
        plt.tight_layout()
