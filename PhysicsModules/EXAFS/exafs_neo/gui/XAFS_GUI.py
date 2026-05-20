"""
Authors    Matthew Adas, Miu Lun(Andy) Lau*, Jeffrey Terry, Min Long
Email      madas@hawk.iit.edu, andylau@u.boisestate.edu, jterry@agni.phys.iit.edu, minlong@boisestate.edu
Version    0.2.2
Date       Sep 6, 2021

Please start the program within the "gui" directory, or "select file" buttons won't start the user in EXAFS.
"Select Directory" button starts the user in EXAFS/path_files/Cu
"""

import ast
import configparser

# Import Basic libraries
import glob
import os
import subprocess

# Import Tkinter
import tkinter as tk

# import queue
from tkinter import (
    ttk,
    Tk,
    N,
    W,
    E,
    S,
    StringVar,
    IntVar,
    DoubleVar,
    BooleanVar,
    filedialog,
    messagebox,
    PhotoImage,
    LabelFrame,
    scrolledtext,
)
from tkinter.font import Font

import numpy as np
from psutil import cpu_count

from PhysicsModules.EXAFS.exafs_neo.helper import Bcolors
from PhysicsModules.EXAFS.exafs_neo.utils_mapping import (
    neomutator_str2int,
    neoselector_str2int,
    neocrossover_str2int,
    neosolver_str2int,
)

os.environ["NUMEXPR_MAX_THREADS"] = str(cpu_count())
# import Larch
import re
import larch

# import matplotlib for figures
import matplotlib

matplotlib.use("TkAgg")

# Import custom libraries
from PhysicsModules.EXAFS.exafs_neo.gui.Background_plot import BKG_plot
from PhysicsModules.EXAFS.exafs_neo.gui.Analysis_plot import Analysis_Plot
from PhysicsModules.EXAFS.exafs_neo.gui.feff_folder_larch import *
from PhysicsModules.EXAFS.exafs_neo._version import __version__
from PhysicsModules.EXAFS.exafs_neo.gui.Misc_Function import natural_keys, get_platform, check_sabcor_folder


class App:
    """
    Start of the applications
    """

    def __init__(self):

        self.__version__ = "0.2.1"
        # Larch
        self.mylarch = larch.Interpreter()
        self.root = Tk(className="EXAFS Neo GUI")
        self.root.wm_title("Graphical User Interface for EXAFS Analysis (Beta)")
        # Standard default geometry
        self.root.geometry("750x500")
        self.padx = 5
        self.pady = 3

        self.os = get_platform()
        base_folder = os.path.dirname(os.path.join(os.getcwd(), __file__))
        icon_loc = os.path.join(base_folder, "media/icon.png")

        img = PhotoImage(file=icon_loc)
        self.root.tk.call("wm", "iconphoto", self.root._w, img)

        # Set default font
        self.entryFont = Font(family="TkFixedFont", size=10)
        self.menuFont = Font(family="TkMenuFont", size=10)
        self.labelFont = Font(family="TkTextFont", size=11)
        # Generate the mainframe
        self.mainframe = ttk.Notebook(self.root, height=40, padding="5")
        self.mainframe.grid(column=0, row=0, sticky=(N, W, E, S))

        # Initialize variable
        self.initialize_var()
        self.initialize_tab()
        self.Build_tabs()

        # Initialize sabcor
        check_sabcor_folder(self.os)

    def initialize_var(self) -> None:
        """
        Initalize all possible variables in the gui.

        Future additions starts here.
        """
        # Inputs
        self.data_file = StringVar(self.root, "Please choose file/s")
        self.temp_data_file = StringVar(self.root, "Please choose file/s")
        self.output_file = StringVar(self.root, "Please choose a file")
        self.ncomp = IntVar(self.root, 0)  # Number of compounds
        self.feff_file = StringVar(self.root, "Please choose a directory")
        self.series = BooleanVar(self.root, False)

        # Populations
        self.populations = IntVar(self.root, 1000)
        self.num_gen = IntVar(self.root, 100)
        self.best_sample = IntVar(self.root, 20)
        self.lucky_few = IntVar(self.root, 20)

        # Mutations
        self.chance_of_mutation = IntVar(self.root, 20)
        self.orginal_chance_of_mutation = IntVar(self.root, 20)
        self.chance_of_mutation_e0 = IntVar(self.root, 20)
        self.mutation_options = StringVar(self.root, "Mutate Per Individual")

        # Selection
        self.selection_options = StringVar(self.root, "Roulette Wheel")

        # Crossover
        self.crossover_options = StringVar(self.root, "Uniform Crossover")

        # Solver
        self.solver_options = StringVar(self.root, "Genetic Algorithm with Rechenberg")

        # Paths
        self.individual_path = BooleanVar(self.root, True)
        self.path_range = IntVar(self.root, 20)
        temp_list = ",".join(np.char.mod("%i", np.arange(1, self.path_range.get() + 1)))
        self.path_list = StringVar(self.root, temp_list)
        self.path_optimize = BooleanVar(self.root, False)
        self.path_optimize_pert = DoubleVar(self.root, 0.01)
        self.path_optimize_only = BooleanVar(self.root, False)

        # Larch Paths
        self.kmin = DoubleVar(self.root, 2.5)
        self.kmax = DoubleVar(self.root, 15)
        self.k_weight = DoubleVar(self.root, 2)
        self.delta_k = DoubleVar(self.root, 0.05)
        self.r_bkg = DoubleVar(self.root, 1.0)
        self.bkg_kw = DoubleVar(self.root, 1.0)
        self.bkg_kmax = DoubleVar(self.root, 15)

        # Outputs
        self.print_graph = BooleanVar(self.root, False)
        self.num_output_paths = BooleanVar(self.root, True)
        self.steady_state_exit = BooleanVar(self.root, True)

        # Pertubutuions
        self.n_ini = IntVar(self.root, 100)

        # Sabcor initalizes variables
        self.sabcor_input_file = StringVar(self.root, "Please choose a file")
        self.sabcor_executable = StringVar(
            self.root, "Please choose the sabcor executable"
        )
        self.sabcor_toggle = BooleanVar(self.root, False)

    def initialize_tab(self) -> None:
        """
        Initialize tab for the main frame. more tabs
        """

        s = ttk.Style()
        s.configure("TNotebook.Tab", font=("TkHeadingFont", "11"))

        height = 1
        self.tab_Inputs = tk.Frame(self.mainframe, height=height)
        self.tab_Populations = tk.Frame(self.mainframe, height=height)
        self.tab_Paths = tk.Frame(self.mainframe, height=height)
        self.tab_Mutation = tk.Frame(self.mainframe, height=height)
        self.tab_Larch = tk.Frame(self.mainframe, height=height)
        self.tab_Bkg = tk.Frame(self.mainframe, height=height)
        self.tab_Output = tk.Frame(self.mainframe, height=height)
        self.tab_Analysis = tk.Frame(self.mainframe, height=height)
        self.tab_Expert = tk.Frame(self.mainframe, height=height)

        # Add all the secondary tabs on the top window
        self.mainframe.add(self.tab_Inputs, text="Inputs")
        self.mainframe.add(self.tab_Populations, text="Populations")
        self.mainframe.add(self.tab_Paths, text="Paths")
        self.mainframe.add(self.tab_Mutation, text="Mutations")
        self.mainframe.add(self.tab_Larch, text="Larch")
        self.mainframe.add(self.tab_Bkg, text="Background Plots")
        self.mainframe.add(self.tab_Output, text="Outputs")
        self.mainframe.add(self.tab_Analysis, text="Analysis")
        self.mainframe.add(self.tab_Expert, text="Expert")
        self.mainframe.grid(
            row=0,
            column=0,
            columnspan=5,
            padx=self.padx,
            pady=self.pady,
            sticky=E + W + N + S,
        )

    def read_ini(self, file_name: Union[str, Path]) -> None:
        config = configparser.ConfigParser()
        config.read(file_name)

        self.ncomp.set(config.getint("Inputs", "num_compounds", fallback=1))
        self.data_file.set(config.get("Inputs", "csv_file"))
        self.temp_data_file.set(config.get("Inputs", "csv_file"))

        self.output_file.set(config.get("Inputs", "output_file"))
        # FEFF file(s) may be a single entry or comma-separated list for
        # multiple components. Update both the aggregated `self.feff_file`
        # and individual per-component variables (created by
        # `gen_feff_folder`) so the GUI entries reflect the ini.
        feff_val = config.get(
            "Inputs", "feff_file", fallback="Please choose a directory"
        )
        self.feff_file.set(feff_val)
        # If the per-component feff variables exist, populate them as well.
        if hasattr(self, "feff_file_list") and self.feff_file_list:
            # split on commas and strip whitespace
            parts = [p.strip() for p in feff_val.split(",") if p.strip()]
            for i in range(min(len(parts), len(self.feff_file_list))):
                try:
                    self.feff_file_list[i].set(parts[i])
                except Exception:
                    # defensive: ignore if something unexpected is present
                    pass

        self.series.set(config.getboolean("Inputs", "csv_series", fallback=False))

        self.populations.set(config.getint("Populations", "population"))
        self.num_gen.set(config.getint("Populations", "num_gen"))
        self.best_sample.set(config.getint("Populations", "best_sample"))
        self.lucky_few.set(config.getint("Populations", "lucky_few"))

        self.chance_of_mutation.set(config.getint("Mutations", "chance_of_mutation"))
        self.orginal_chance_of_mutation.set(
            config.getint("Mutations", "original_chance_of_mutation")
        )
        self.chance_of_mutation_e0.set(
            config.getint("Mutations", "chance_of_mutation_e0")
        )
        self.mutation_options.set(
            config.get("Mutations", "mutated_options", fallback=0)
        )
        # TODO: this should be in own separate dictionary
        self.selection_options.set(
            config.get("Mutations", "selection_options", fallback=0)
        )

        self.crossover_options.set(
            config.get("Crossover", "crossover_options", fallback=0)
        )

        self.solver_options.set(config.get("Solver", "solver_options", fallback=1))

        self.individual_path.set(config.getboolean("Paths", "individual_path"))
        self.path_range.set(config.getint("Paths", "path_range"))
        self.path_list.set(config.get("Paths", "path_list"))

        if hasattr(self, "path_file_list") and self.path_file_list:
            input_string = self.path_list.get()
            if len(self.path_file_list) > 1:
                groups = re.findall(r"\[([^\]]*)\]", input_string)
                if not groups:
                    # no brackets: treat entire string as single component
                    groups = [input_string] if input_string else []
                groups = [g.strip() for g in groups]
                for i in range(len(groups)):
                    self.path_file_list[i].set(groups[i])
            else:
                self.path_file_list[0].set(input_string)

        self.path_optimize.set(
            config.getboolean("Paths", "path_optimize", fallback=False)
        )

        self.kmin.set(config.getfloat("Larch_Paths", "kmin"))
        self.kmax.set(config.getfloat("Larch_Paths", "kmax"))
        self.k_weight.set(config.getfloat("Larch_Paths", "kweight"))
        self.delta_k.set(config.getfloat("Larch_Paths", "deltak"))
        self.r_bkg.set(config.getfloat("Larch_Paths", "rbkg"))
        self.bkg_kw.set(config.getfloat("Larch_Paths", "bkgkw"))
        self.bkg_kmax.set(config.getfloat("Larch_Paths", "bkgkmax"))

        self.print_graph.set(config.getboolean("Outputs", "print_graph"))
        self.num_output_paths.set(config.getboolean("Outputs", "num_output_paths"))
        self.steady_state_exit.set(config.getboolean("Outputs", "steady_state_exit"))

    def description_tabs(self, arr, tabs, sticky=(W, E), row=None) -> None:
        # Rows = index of rows
        description_list = []
        if row is not None:
            assert len(row) == len(arr)
        for i, inputs in enumerate(arr):
            entry = ttk.Label(tabs, text=inputs, font=self.labelFont)
            if row is not None:
                k = row[i]
            else:
                k = i
            entry.grid_configure(
                column=0, row=k, sticky=sticky, padx=self.padx, pady=self.pady
            )
            description_list.append(entry)
        if description_list:
            return description_list
        else:
            return []

    def Build_tabs(self) -> None:
        """
        Build tabs
        """
        # Building all arrays
        self.Build_global()
        self.Build_inputs_tab()
        self.Build_population_tab()
        self.Build_path_tab()
        self.Build_mutations_tab()
        self.Build_larch_tab()
        self.Build_background_tab()
        self.Build_output_tab()
        self.Build_analysis_tab()
        self.Build_expert_tab()

        self.mainframe.grid_rowconfigure(0, weight=1)
        self.mainframe.grid_columnconfigure(0, weight=1)
        self.mainframe.grid_columnconfigure(1, weight=1)

    def write_ini(self, filename: Union[str, Path]) -> None:
        """
        Write the ini for the specific file
        """
        inputs = "[Inputs] \nnum_compounds = {compound} \ncsv_file = {data} \noutput_file = {out} \nfeff_file = {feff}\ncsv_series = {series}".format(
            compound=str(self.ncomp.get()),
            data=str(self.data_file.get()),
            out=str(self.output_file.get()),
            feff=str(self.feff_file.get()),
            series=str(self.series.get()),
        )

        populations = "\n\n[Populations] \npopulation = {pop} \nnum_gen = {numgen} \nbest_sample = {best} \nlucky_few = {luck}".format(
            pop=str(self.populations.get()),
            numgen=str(self.num_gen.get()),
            best=str(self.best_sample.get()),
            luck=str(self.lucky_few.get()),
        )

        mutations = "\n\n[Mutations] \nchance_of_mutation = {chance} \noriginal_chance_of_mutation = {original} \nchance_of_mutation_e0 = {e0} \nmutated_options = {opt}".format(
            chance=str(self.chance_of_mutation.get()),
            original=str(self.orginal_chance_of_mutation.get()),
            e0=str(self.chance_of_mutation_e0.get()),
            opt=str(neomutator_str2int(self.mutation_options.get())),
        )

        selector = "\n\n[Selector] \nselection_options = {sel}".format(
            sel=str(neoselector_str2int(self.selection_options.get()))
        )

        crossover = "\n\n[Crossover] \ncrossover_options = {croOpt}".format(
            croOpt=str(neocrossover_str2int(self.crossover_options.get()))
        )

        solver = "\n\n[Solver] \nsolver_options = {solOpt}".format(
            solOpt=str(neosolver_str2int(self.solver_options.get()))
        )

        paths = "\n\n[Paths] \nindividual_path = {tf}  \npath_range = {range} \npath_list = {list} \npath_optimize = {optimize} \npath_optimize_percent = {optimize_pert} \noptimize_only = {optimize_only}".format(
            tf=str(self.individual_path.get()),
            range=str(self.path_range.get()),
            list=str(self.path_list.get().replace(" ", "")),
            optimize=str(self.path_optimize.get()),
            optimize_pert=str(self.path_optimize_pert.get()),
            optimize_only=str(self.path_optimize_only.get()),
        )

        larch_paths = "\n\n[Larch_Paths] \nkmin = {min} \nkmax = {max} \nkweight = {weight} \ndeltak = {delk} \nrbkg = {rb} \nbkgkw = {bk} \nbkgkmax = {bmax}".format(
            min=self.kmin.get(),
            max=self.kmax.get(),
            weight=self.k_weight.get(),
            delk=self.delta_k.get(),
            rb=self.r_bkg.get(),
            bk=self.bkg_kw.get(),
            bmax=self.bkg_kmax.get(),
        )

        outputs = "\n\n[Outputs] \nprint_graph = {pg} \nnum_output_paths = {outpath}\nsteady_state_exit = {steady}".format(
            pg=self.print_graph.get(),
            outpath=self.num_output_paths.get(),
            steady=self.steady_state_exit.get(),
        )

        with open(filename, "w") as writer:
            writer.write(str(inputs))
            writer.write(str(populations))
            writer.write(str(solver))
            writer.write(str(mutations))
            writer.write(str(crossover))
            writer.write(str(selector))
            writer.write(str(paths))
            writer.write(str(larch_paths))
            writer.write(str(outputs))

    def generate_ini(self) -> None:
        ini_file = filedialog.asksaveasfilename(
            initialdir=os.getcwd(),
            title="Choose output ini file",
            filetypes=[("ini files", "*.ini")],
        )
        if ini_file is None:
            return
        if not isinstance(ini_file, tuple):
            if len(ini_file) != 0:
                self.write_ini(ini_file)
                messagebox.showinfo(
                    "", "Ini file written to {fileloc}".format(fileloc=ini_file)
                )

    def stop_term(self) -> None:
        if hasattr(self, "proc"):
            self.proc.kill()

    def run_term(self, file: Union[str, Path] = "test_temp.i") -> None:
        """
        if hasattr(self,'terminal'):
            # Need to close previous threads
            self.terminal.destroy()
        # command = ['exafs_neo','-i','test.ini']
        self.Write_ini('test_temp.i')
        command = 'exafs_neo -i test_temp.i'.split(' ')
        # print(command)
        self.terminal = Console(self.txtbox,command)
        """
        self.stop_term()

        # command = 'exafs_neo -i test_temp.i'
        if self.sabcor_toggle.get():
            # excut_path = check_executable()
            import sabcor

            execut_path = self.sabcor_executable.get()
            sabcor.check_executable(paths=execut_path)
            params = sabcor.read_sab(self.sabcor_input_file.get())
            sabcor.write_sab(params)

            sabcor.call_executable(execut_path, self.data_file.get())
            sabcor.edited_final_header(self.data_file.get())
            post_sabcor_file = os.path.splitext(self.data_file.get())
            self.data_file.set(post_sabcor_file[0] + "_sac" + post_sabcor_file[1])

        self.write_ini(str(Path.cwd() / "test_temp.i"))

        command = ["exafs_neo", "-i", file]
        if self.os == "Windows":
            print(" ".join(command))
            self.proc = subprocess.Popen(" ".join(command), shell=True)
        else:
            self.proc = subprocess.Popen("exec " + " ".join(command), shell=True)

    def run_ini(self, file="test_temp.i") -> None:

        command = ["exafs_neo", "-i", file]
        if self.os == "Windows":
            self.proc = subprocess.Popen("call " + " ".join(command), shell=True)
        else:
            self.proc = subprocess.Popen("exec " + " ".join(command), shell=True)
        self.proc.wait()

    def Build_global(self):
        """
        Create global generate ini
        """

        def about_citation() -> None:
            """
            Create about popup
            """
            popup = tk.Toplevel()
            popup.wm_title("About: Ver: " + str(self.__version__))

            popup.resizable(False, False)
            cite = tk.Label(popup, text="Citation:", font="TkTextFont")
            cite.grid(column=0, row=0, sticky=N, padx=self.padx, pady=self.pady)
            citation = scrolledtext.ScrolledText(
                popup, width=75, height=10, font="TkTextFont"
            )
            citation.grid(column=0, row=1, sticky=N, padx=self.padx, pady=self.pady)
            file_path = os.path.dirname(__file__)
            with open(os.path.join(file_path, "media/Citation")) as f:
                citation.insert(tk.END, f.read())
            citation.configure(state="disabled")

            license_label = tk.Label(popup, text="License:", font="TkTextFont")
            license_label.grid(
                column=0, row=2, sticky=N, padx=self.padx, pady=self.pady
            )
            license_text = scrolledtext.ScrolledText(popup, width=75, font="TkTextFont")
            license_text.grid(column=0, row=3, sticky=N, padx=self.padx, pady=self.pady)
            with open(os.path.join(file_path, "media/LICENSE")) as f:
                license_text.insert(tk.END, f.read())
            license_text.configure(state="disabled")

            B1 = ttk.Button(popup, text="Close", command=popup.destroy)
            B1.grid(column=0, row=4, padx=self.padx, pady=self.pady)

            popup.grid_columnconfigure((0, 2), weight=1)
            popup.grid_rowconfigure((0, 2), weight=1)

            popup.protocol("WM_DELETE_WINDOW", popup.destroy)

        def load_ini() -> None:
            file_name = filedialog.askopenfilename(
                initialdir=os.getcwd(),
                title="Choose ini file",
                filetypes=(("ini files", "*.ini"),),
            )
            if not file_name:
                # Raise Warning
                self.ini_file.set("Please choose file/s")
                self.temp_ini_file.set("Please choose file/s")
            else:
                self.read_ini(file_name)

        self.generate_button = tk.Button(
            self.root, text="Generate Input", command=self.generate_ini
        )
        self.generate_button.grid(
            column=4, row=2, sticky=E, padx=self.padx, pady=self.pady
        )

        self.run_button = tk.Button(self.root, text="Run", command=self.run_term)
        self.run_button.grid(
            column=2, row=2, columnspan=1, sticky=E, padx=self.padx, pady=self.pady
        )

        self.stop_button = tk.Button(self.root, text="Stop", command=self.stop_term)
        self.stop_button.grid(
            column=3, row=2, columnspan=1, sticky=W, padx=self.padx, pady=self.pady
        )

        if self.os == "Windows":
            self.stop_button.config(state="disabled")

        self.about_button = tk.Button(self.root, text="?", command=about_citation)
        self.about_button.grid(
            column=0, row=2, columnspan=1, sticky=W, padx=self.padx, pady=self.pady
        )

        self.load_ini_button = tk.Button(self.root, text="Load Ini", command=load_ini)
        self.load_ini_button.grid(
            column=1, row=2, columnspan=1, sticky=W, padx=self.padx, pady=self.pady
        )

        self.root.grid_columnconfigure(2, weight=1)
        self.root.grid_columnconfigure(3, weight=1)

        self.root.grid_rowconfigure(0, weight=1)

        # Create an empty frame
        self.label_frame = LabelFrame(self.root, text="Terminal", padx=5, pady=5)
        self.label_frame.grid(
            column=0,
            row=1,
            columnspan=5,
            padx=self.padx,
            pady=self.pady,
            sticky=E + W + N + S,
        )

        # Create the textbox
        self.label_frame.rowconfigure(0, weight=1)
        self.label_frame.columnconfigure(0, weight=1)
        self.txtbox = scrolledtext.ScrolledText(self.label_frame, width=40, height=10)
        self.txtbox.grid(row=0, column=0, sticky=E + W + N + S)

    def Build_inputs_tab(self) -> None:
        arr_input = ["Data file", "Output File", "Number of FEFF folder"]
        # FEFF Folder
        self.description_tabs(arr_input, self.tab_Inputs, row=[0, 1, 3])

        self.tab_Inputs.grid_columnconfigure(1, weight=1)

        entry_data_file = ttk.Combobox(
            self.tab_Inputs, textvariable=self.temp_data_file, font=self.entryFont
        )
        entry_data_file.grid(
            column=1, row=0, sticky=(W, E), padx=self.padx, pady=self.pady
        )

        entry_output_file = tk.Entry(
            self.tab_Inputs, textvariable=self.output_file, font=self.entryFont
        )
        entry_output_file.grid(
            column=1, row=1, sticky=(W, E), padx=self.padx, pady=self.pady
        )

        separator = ttk.Separator(self.tab_Inputs, orient="horizontal")
        separator.grid(column=0, row=2, columnspan=4, sticky=W + E, padx=self.padx)

        comp_list = list(range(1, 6))
        entry_ncomp = ttk.Combobox(
            self.tab_Inputs,
            width=7,
            values=comp_list,
            textvariable=self.ncomp,
            font=self.entryFont,
        )
        entry_ncomp.grid(column=1, row=3, sticky=(W, E), padx=self.padx)

        def select_data_file():

            file_name = filedialog.askopenfilenames(
                initialdir=os.getcwd(),
                title="Choose xmu/csv",
                filetypes=(
                    ("xmu files", "*.xmu"),
                    ("csv files", "*.csv"),
                    ("all files", "*.*"),
                ),
            )
            if not file_name:
                self.data_file.set("Please choose file/s")
                self.temp_data_file.set("Please choose file/s")
            else:
                if isinstance(file_name, tuple):
                    if len(file_name) == 1:
                        self.data_file.set(str(file_name[0]))
                        self.temp_data_file.set(str(file_name[0]))
                    else:
                        file_list = list(self.root.splitlist(file_name))
                        self.data_file.set(" ".join(file_list))
                        entry_data_file["value"] = file_list
                        self.series.set(True)

        def select_output_file() -> None:
            file_name = filedialog.asksaveasfilename(
                initialdir=os.getcwd(),
                title="Choose data_file",
                filetypes=(("csv files", "*.csv"), ("all files", "*.*")),
            )
            if not file_name:
                self.output_file.set("Please choose a file")
            else:
                self.output_file.set(file_name)

        def feff_trace(var, indx, mode):
            # Todo: Need to add assertion to make sure every folder is not normal
            raw_str = []
            for i in range(len(self.feff_file_list)):
                raw_str.append(self.feff_file_list[i].get())
            raw_feff = ",".join(raw_str)
            self.feff_file.set(raw_feff)

        def gen_feff_folder(var, indx, mode):
            ncomp = self.ncomp.get()
            self.feff_file_toggle = [True for i in range(self.ncomp.get())]
            k = 4
            try:
                self.arr_feff
                # self.input_list
            except AttributeError:
                pass
            else:
                for i in range(len(self.feff_input_list)):
                    self.feff_input_list[i].destroy()
                    self.feff_button_list[i].destroy()
                    self.feff_description_list[i].destroy()
            self.feff_file_list = []
            self.arr_feff = []
            arr_row = []
            # parameter sets
            self.feff_input_list = []
            self.feff_button_list = []
            self.feff_description_list = []
            # Need to remove previous one
            for i in range(ncomp):
                self.arr_feff.append("FEFF Folder (" + str(i + 1) + ")")
                arr_row.append(i + k)

                # Generate a list of variables. If `self.feff_file` already
                # contains FEFF entries from a loaded ini, initialize each
                # per-component variable from that value so we don't overwrite
                # the user's selection when the GUI regenerates these fields.
                temp_var = StringVar(self.root, "Please choose a directory")
                # If aggregated feff_file has been set (comma-separated for
                # multiple components), use the corresponding part.
                try:
                    feff_agg = self.feff_file.get()
                except Exception:
                    feff_agg = ""
                if feff_agg and feff_agg != "Please choose a directory":
                    parts = [p.strip() for p in feff_agg.split(",") if p.strip()]
                    if i < len(parts):
                        temp_var.set(parts[i])

                self.feff_file_list.append(temp_var)

                # Add trace back
                temp_var.trace_add("write", feff_trace)
                # Setup each entry
                temp_entry = tk.Entry(
                    self.tab_Inputs, textvariable=temp_var, font=self.entryFont
                )

                temp_entry.grid(
                    column=1, row=k + i, sticky=(W, E), padx=self.padx, pady=self.pady
                )
                self.feff_input_list.append(temp_entry)

                # setup each button
                temp_button = ttk.Button(
                    self.tab_Inputs,
                    text="Choose",
                    command=lambda x=i: select_feff_folder(self.feff_file_list[x]),
                    style="my.TButton",
                )
                temp_button.grid(
                    column=3, row=k + i, sticky=W, padx=self.padx, pady=self.pady
                )
                self.feff_button_list.append(temp_button)

            self.feff_description_list = self.description_tabs(
                self.arr_feff, self.tab_Inputs, row=arr_row
            )

        def check_feff_folder(folder: Union[str, Path]) -> None:
            feff_inp_exists = False
            num_feff_file = 0
            folder_Path = Path(folder)
            for i in folder_Path.iterdir():
                file = Path(i).name
                if fnmatch.fnmatch(file, "feff????.dat"):
                    num_feff_file += 1

            if num_feff_file == 0:
                # Check no scattering paths
                for i in folder_Path.iterdir():
                    file = Path(i).name
                    if fnmatch.fnmatch(file, "feff.inp"):
                        feff_inp_exists = True
                        feff_inp_loc = file
                abs_feff_inp_loc = folder_Path.joinpath("feff.inp")
                create_feff_folder(abs_feff_inp_loc)
            else:
                self.txtbox.insert(tk.END, "No feff.inp in folder\n")

        def select_feff_folder(var) -> None:
            # select folder name
            folder_name = filedialog.askdirectory(
                initialdir=os.getcwd(), title="Select folder"
            )

            if not folder_name:
                var.set("Please choose a directory")
            else:
                check_feff_folder(folder_name)
                folder_name = os.path.join(folder_name, "feff")
                var.set(folder_name)

        # add trace back
        self.ncomp.trace_add("write", gen_feff_folder)
        self.ncomp.set(1)

        button_data_file = ttk.Button(
            self.tab_Inputs, text="Choose", command=select_data_file, style="my.TButton"
        )
        button_data_file.grid(column=3, row=0, sticky=W, padx=self.padx, pady=self.pady)

        button_output_file = ttk.Button(
            self.tab_Inputs,
            text="Choose",
            command=select_output_file,
            style="my.TButton",
        )
        button_output_file.grid(
            column=3, row=1, sticky=W, padx=self.padx, pady=self.pady
        )

    def Build_population_tab(self):
        arr_pop = [
            "Population",
            "Number of Generations",
            "Best Individuals(%)",
            "Lucky Survivor(%)",
        ]
        self.description_tabs(arr_pop, self.tab_Populations)

        self.tab_Populations.grid_columnconfigure(1, weight=1)

        entry_population = tk.Entry(
            self.tab_Populations, textvariable=self.populations, font=self.entryFont
        )
        entry_population.grid(column=1, row=0, sticky=(W, E), padx=self.padx)

        entry_num_gen = tk.Entry(
            self.tab_Populations, textvariable=self.num_gen, font=self.entryFont
        )
        entry_num_gen.grid(column=1, row=1, sticky=(W, E), padx=self.padx)

        entry_best_sample = tk.Entry(
            self.tab_Populations, textvariable=self.best_sample, font=self.entryFont
        )
        entry_best_sample.grid(column=1, row=2, sticky=(W, E), padx=self.padx)

        entry_lucky_few = tk.Entry(
            self.tab_Populations, textvariable=self.lucky_few, font=self.entryFont
        )
        entry_lucky_few.grid(column=1, row=3, sticky=(W, E), padx=self.padx)

    def Build_path_tab(self):
        """
        Build path tabs
        """

        arr_paths = ["Path Optimize", "Multi-Path Optimize", "Path Optimize Percentage"]
        self.description_tabs(arr_paths, self.tab_Paths, row=[0, 1, 2])
        """
        TODO: Disable this checkbox calc for individual paths
        def checkbox_individual_paths():
            if self.individual_path.get() == True:
                pass
                entry_path_range.config(state="disabled")
                entry_path_list.config(state="normal")

            else:
                pass
                entry_path_range.config(state="normal")
                entry_path_list.config(state="disabled")
        """

        # def path_list_cb(var, indx, mode):
        #     """
        #     Path_list call-back for the number of lists
        #     """
        #     if self.individual_path.get()== True:
        #         self.path_list.set(self.path_list_call.get())
        #         counts = 0
        #         test_list = self.path_list.get().split(",")
        #         for i in range(len(test_list)):
        #             if test_list[i] != '':
        #                 counts += 1
        #         self.path_range_call.set(counts)

        # entry_individual_paths = ttk.Checkbutton(self.tab_Paths,
        #     variable = self.individual_path,command = checkbox_individual_paths)
        # entry_individual_paths.grid(column=1,row=0,sticky=(W),padx=self.padx)

        # self.path_range_call = StringVar(self.tab_Paths,str(self.path_range.get()))
        # # self.path_range_call.trace_add("write",path_range_cb)
        # entry_path_range = ttk.Entry(self.tab_Paths, width=7,
        #         textvariable=self.path_range_call, font=self.entryFont)
        # entry_path_range.grid(column=1, row=1, sticky=(W),padx=self.padx)

        # self.path_list_call = StringVar(self.tab_Paths,self.path_list.get())
        # self.path_list_call.trace_add('write',path_list_cb)
        # entry_path_list = ttk.Entry(self.tab_Paths,
        #         textvariable=self.path_list_call, font=self.entryFont)
        # entry_path_list.grid(column=1, row=2, sticky=(W, E),padx=self.padx)
        # entry_path_list.config(state='disabled')

        # def path_range_cb(var, indx, mode):
        #     if self.individual_path.get()==False:
        #         if self.path_range_call.get() == '':
        #             int_test = 0
        #         else:
        #             int_test = int(float(self.path_range_call.get()))
        #         self.path_range_call.set(int_test)
        #         self.path_range.set(int_test)
        #
        #         custom_path_list = ",".join(np.char.mod('%i', np.arange(1,self.path_range.get()+1)))
        #         self.path_list_call.set(custom_path_list)
        #         self.path_list.set(custom_path_list)
        # def test_str(var):

        def output_var_trace(var, indx, mode):
            """
            Updates the variables whenever it gets changes
            """
            raw_str = []
            # compounds
            if self.ncomp.get() > 1:
                for i in range(len(self.path_file_list)):
                    var = self.path_file_list[i].get()
                    if var == str(0) or var == "":
                        base_str = ""
                        # del the corresponding one from the final output
                        self.feff_file_toggle[i] = False
                    else:
                        base_str = "[" + var + "]"
                        self.feff_file_toggle[i] = True
                    raw_str.append(base_str)

                # call back for path and feff
                temp_path_list = []
                temp_feff_list = []
                for i in range(len(self.path_file_list)):
                    if self.feff_file_toggle[i] == True:
                        temp_path_list.append(raw_str[i])
                        temp_feff_list.append(self.feff_file_list[i].get())
                raw_path = ",".join(temp_path_list)
                raw_feff = ",".join(temp_feff_list)
                self.path_list.set(raw_path)
                self.feff_file.set(raw_feff)
                self.individual_path.set(True)
            # else calculate compounds
            else:
                self.path_list.set(self.path_file_list[0].get())

            # print("-----------------------")
            # print(self.feff_file_toggle)
            # print(self.feff_file.get())
            # print(self.path_list.get())

        def ncomp_trace_cb(var, indx, mode):
            """
            Trace the number of components, if changes, recomputed
            """
            ncomp = self.ncomp.get()
            k = 5  # Starting components

            # Check if arr_pathlist exists
            try:
                self.arr_pathlist
            except AttributeError:
                pass
            else:
                # Destory all in the current list
                for i in range(len(self.path_input_list)):
                    self.path_description_list[i].destroy()
                    self.path_input_list[i].destroy()
            # Create the list
            self.arr_pathlist = []
            self.path_input_list = []
            self.path_file_list = []
            arr_row = []
            for i in range(ncomp):
                self.arr_pathlist.append("Pathlist Folder (" + str(i + 1) + ")")
                arr_row.append(i + k)
                # Create the entry variable
                temp_path_var = StringVar(self.root, "Path_list")
                temp_path_var.trace_add("write", output_var_trace)
                # temp_path_var.trace_add('write',)
                self.path_file_list.append(temp_path_var)
                # Create the entry
                temp_entry = tk.Entry(
                    self.tab_Paths,
                    textvariable=self.path_file_list[i],
                    font=self.entryFont,
                )
                # lock to the grid
                temp_entry.grid(
                    column=1,
                    columnspan=2,
                    row=k + i,
                    sticky=(W, E),
                    padx=self.padx,
                    pady=self.pady,
                )

                self.path_input_list.append(temp_entry)

            # create the descriptions for it.
            self.path_description_list = self.description_tabs(
                self.arr_pathlist, self.tab_Paths, row=arr_row
            )

        entry_path_optimize = ttk.Checkbutton(
            self.tab_Paths, variable=self.path_optimize
        )
        entry_path_optimize.grid(column=1, row=0, sticky=(W, E), padx=self.padx)

        entry_path_optimize_only = ttk.Checkbutton(
            self.tab_Paths, variable=self.path_optimize_only
        )
        entry_path_optimize_only.grid(column=1, row=1, sticky=(W, E), padx=self.padx)

        entry_optimize_perc = ttk.Entry(
            self.tab_Paths, textvariable=self.path_optimize_pert, font=self.entryFont
        )
        entry_optimize_perc.grid(column=1, row=2, sticky=(W, E), padx=self.padx)

        separator = ttk.Separator(self.tab_Paths, orient="horizontal")
        separator.grid(column=0, row=3, columnspan=4, sticky=(W, E), padx=self.padx)
        self.tab_Paths.columnconfigure(4, weight=1)

        # add the call back for it to modify if it get changes
        self.ncomp.trace_add("write", ncomp_trace_cb)

        # inital default value
        self.ncomp.set(1)

    def Build_mutations_tab(self):
        arr_muts = [
            "Mutation Chance (%)",
            "Original chance of mutation (%)",
            "E0 Mutation Chance (%)",
            "Mutation Options",
        ]
        self.description_tabs(arr_muts, self.tab_Mutation)

        mut_list = list(range(101))
        # self.tab_Inputs.grid_columnconfigure(1,weight=1)

        entry_chance_of_mutation = ttk.Combobox(
            self.tab_Mutation,
            width=7,
            values=mut_list,
            textvariable=self.chance_of_mutation,
            font=self.entryFont,
        )
        entry_chance_of_mutation.grid(column=1, row=0, sticky=(W, E), padx=self.padx)

        entry_chance_of_mutation = ttk.Combobox(
            self.tab_Mutation,
            width=7,
            values=mut_list,
            textvariable=self.orginal_chance_of_mutation,
            font=self.entryFont,
        )
        entry_chance_of_mutation.grid(column=1, row=1, sticky=(W, E), padx=self.padx)

        entry_chance_of_mutation = ttk.Combobox(
            self.tab_Mutation,
            width=7,
            values=mut_list,
            textvariable=self.chance_of_mutation_e0,
            font=self.entryFont,
        )
        entry_chance_of_mutation.grid(column=1, row=2, sticky=(W, E), padx=self.padx)

        mut_list = [
            "Mutate Per Individual",
            "Mutate Per Path",
            "Mutate Per Trait",
            "Mutate Metropolis",
            "Mutate Bounded Per Range",
            "Mutate Differential Evolution",
        ]
        entry_chance_of_mutation = ttk.Combobox(
            self.tab_Mutation,
            width=7,
            values=mut_list,
            textvariable=self.mutation_options,
            font=self.entryFont,
        )
        entry_chance_of_mutation.grid(column=1, row=3, sticky=(W, E), padx=self.padx)

        self.tab_Mutation.columnconfigure(1, weight=1)

    def Build_larch_tab(self):
        arr_larch = [
            "Kmin",
            "Kmax",
            "Kweight",
            "Delta k",
            "R Bkg",
            "Bkg Kw",
            "Bkg Kmax",
        ]
        self.description_tabs(arr_larch, self.tab_Larch)

        entry_kmin = ttk.Entry(
            self.tab_Larch, textvariable=self.kmin, font=self.entryFont
        )
        entry_kmin.grid(column=1, row=0, sticky=(W, E), padx=self.padx)

        entry_kmax = ttk.Entry(
            self.tab_Larch, textvariable=self.kmax, font=self.entryFont
        )
        entry_kmax.grid(column=1, row=1, sticky=(W, E), padx=self.padx)

        entry_k_weight = ttk.Entry(
            self.tab_Larch, textvariable=self.k_weight, font=self.entryFont
        )
        entry_k_weight.grid(column=1, row=2, sticky=(W, E), padx=self.padx)

        entry_delta_k = ttk.Entry(
            self.tab_Larch, textvariable=self.delta_k, font=self.entryFont
        )
        entry_delta_k.grid(column=1, row=3, sticky=(W, E), padx=self.padx)

        entry_r_bkg = ttk.Entry(
            self.tab_Larch, textvariable=self.r_bkg, font=self.entryFont
        )
        entry_r_bkg.grid(column=1, row=4, sticky=(W, E), padx=self.padx)

        entry_bkg_kw = ttk.Entry(
            self.tab_Larch, textvariable=self.bkg_kw, font=self.entryFont
        )
        entry_bkg_kw.grid(column=1, row=5, sticky=(W, E), padx=self.padx)

        entry_bkg_kmax = ttk.Entry(
            self.tab_Larch, textvariable=self.bkg_kmax, font=self.entryFont
        )
        entry_bkg_kmax.grid(column=1, row=6, sticky=(W, E), padx=self.padx)

    def Build_background_tab(self):

        def button_bkg_draw():
            # TODO:
            # Makes it so
            if self.sabcor_toggle.get() == True:
                base_name = os.path.basename(self.temp_data_file.get())
                base_filename = base_name.split(".")
                # Check if sac is in:
                if "_sac" not in base_filename[0]:
                    import sabcor

                    execut_path = Path.cwd().parent / "contrib/sabcor/bin/sabcor"
                    sabcor.check_executable(paths=execut_path)
                    params = sabcor.read_sab(self.sabcor_input_file.get())
                    sabcor.write_sab(params)
                    print(self.data_file.get())

                    sabcor.call_executable(execut_path, self.data_file.get())
                    sabcor.edited_final_header(self.data_file.get())
                    post_sabcor_file = os.path.splitext(self.data_file.get())
                    # self.data_file.set(post_sabcor_file[0] + "_sac" + post_sabcor_file[1])
                    self.temp_data_file.set(
                        post_sabcor_file[0] + "_sac" + post_sabcor_file[1]
                    )

                # Check if the temp data file contains any sac_
            bkg_plot.initialize_params(
                self.temp_data_file,
                self.r_bkg,
                self.bkg_kw,
                self.bkg_kmax,
                self.kmin,
                self.kmax,
                self.delta_k,
                self.k_weight,
            )
            bkg_plot.draw_background()

        arr_bkg = ["R Bkg", "K min", "K max", "Bkg Kw", "Bkg Kmax"]
        self.description_tabs(arr_bkg, self.tab_Bkg, sticky=(W, E, N))
        # self.tab_Bkg.grid_rowconfigure(3,weight=1)
        # self.tab_Bkg.grid_rowconfigure(4,weight=1)
        # self.tab_Bkg.grid_rowconfigure(5,weight=1)
        # self.tab_Bkg.grid_rowconfigure(6,weight=1)

        self.tab_Bkg.rowconfigure(8, weight=1)
        self.tab_Bkg.columnconfigure(2, weight=1)

        entry_r_bkg = ttk.Entry(
            self.tab_Bkg, textvariable=self.r_bkg, font=self.entryFont
        )
        entry_r_bkg.grid(
            column=1, row=0, sticky=(W, E, N), padx=self.padx, pady=self.pady
        )

        entry_k_min = ttk.Entry(
            self.tab_Bkg, textvariable=self.kmin, font=self.entryFont
        )
        entry_k_min.grid(
            column=1, row=1, sticky=(W, E, N), padx=self.padx, pady=self.pady
        )

        entry_k_max = ttk.Entry(
            self.tab_Bkg, textvariable=self.kmax, font=self.entryFont
        )
        entry_k_max.grid(
            column=1, row=2, sticky=(W, E, N), padx=self.padx, pady=self.pady
        )

        entry_bkg_kw = ttk.Entry(
            self.tab_Bkg, textvariable=self.bkg_kw, font=self.entryFont
        )
        entry_bkg_kw.grid(
            column=1, row=3, sticky=(W, E, N), padx=self.padx, pady=self.pady
        )

        entry_bkg_kmax = ttk.Entry(
            self.tab_Bkg, textvariable=self.bkg_kmax, font=self.entryFont
        )
        entry_bkg_kmax.grid(
            column=1, row=4, sticky=(W, E, N), padx=self.padx, pady=self.pady
        )

        bkg_plot = BKG_plot(self.tab_Bkg, self.mylarch)

        button_bkg_mu = ttk.Button(
            self.tab_Bkg, text="Plot Background and Mu", command=button_bkg_draw
        )
        button_bkg_mu.grid(
            column=0, row=5, columnspan=2, sticky=W + E, padx=self.padx, pady=self.pady
        )

        button_kspace = ttk.Button(
            self.tab_Bkg, text="K space", command=bkg_plot.draw_kspace
        )
        button_kspace.grid(
            column=0, row=6, columnspan=2, sticky=W + E, padx=self.padx, pady=self.pady
        )

        button_rspace = ttk.Button(
            self.tab_Bkg, text="R space", command=bkg_plot.draw_rspace
        )
        button_rspace.grid(
            column=0, row=7, columnspan=2, sticky=W + E, padx=self.padx, pady=self.pady
        )

    def Build_output_tab(self):

        pop_min = IntVar(self.tab_Output, 100)
        pop_max = IntVar(self.tab_Output, 5001)
        gen_min = IntVar(self.tab_Output, 20)
        gen_max = IntVar(self.tab_Output, 501)
        mut_min = IntVar(self.tab_Output, 20)
        mut_max = IntVar(self.tab_Output, 51)

        pertub_check = IntVar(self.tab_Output, 0)

        def generate_multi_ini():
            # os.chdir("..")
            pop_range = np.arange(pop_min.get(), pop_max.get(), 100)
            gen_range = np.arange(gen_min.get(), gen_max.get(), 5)
            mut_range = np.arange(mut_min.get(), mut_max.get(), 10)
            # exit()
            multi_folder = filedialog.askdirectory(
                initialdir=os.getcwd(), title="Select folder"
            )
            if not multi_folder:
                return
            else:
                og_pop = self.populations.get()
                og_gen = self.num_gen.get()
                og_mut = self.chance_of_mutation.get()
                og_out_file = self.output_file.get()

                for i in range(self.n_ini.get()):
                    pop_select = np.random.choice(pop_range)
                    gen_select = np.random.choice(gen_range)
                    mut_select = np.random.choice(mut_range)
                    base_file = os.path.splitext(og_out_file)[0]
                    self.output_file.set(
                        os.path.join(base_file + "_" + str(i).zfill(3) + ".ini")
                    )
                    self.populations.set(pop_select)
                    self.num_gen.set(gen_select)
                    self.chance_of_mutation.set(mut_select)
                    self.write_ini(multi_folder + "/file_" + str(i).zfill(3) + ".i")
                # set back og ini
                self.populations.set(og_pop)
                self.num_gen.set(og_gen)
                self.chance_of_mutation.set(og_mut)
                self.output_file.set(og_out_file)
            # os.chdir("gui")
            return multi_folder

        def run_multi_ini():
            """
            Run multiple ini file
            """
            folder_loc = generate_multi_ini()

            full_file_list = []
            file_list = glob.glob(folder_loc + "/*.i")
            for i in file_list:
                full_file_list.append(os.path.join(folder_loc, i))

            full_file_list.sort(key=natural_keys)
            for i in full_file_list:
                print(Bcolors.BOLD + str(i) + Bcolors.ENDC)
                self.run_ini(i)

        def checkbox_multi():
            widget_lists = [
                entry_n_ini,
                entry_pertub_pop_min,
                entry_pertub_pop_max,
                entry_pertub_gen_min,
                entry_pertub_gen_max,
                entry_pertub_mut_min,
                entry_pertub_mut_max,
                button_gen_nini,
                button_run_nini,
            ]
            if pertub_check.get() == 0:
                for i in widget_lists:
                    i.config(state="disabled")
            elif pertub_check.get() == 1:
                for i in widget_lists:
                    i.config(state="normal")

        arr_out = ["Print graph", "Steady state exit"]
        self.description_tabs(arr_out, self.tab_Output)

        checkbutton_print_graph = ttk.Checkbutton(self.tab_Output, var=self.print_graph)
        checkbutton_print_graph.grid(column=1, row=0, sticky=W + E, padx=self.padx)

        checkbutton_steady_state = ttk.Checkbutton(
            self.tab_Output, var=self.steady_state_exit
        )
        checkbutton_steady_state.grid(column=1, row=1, sticky=W + E, padx=self.padx)

        # Create separators
        separator = ttk.Separator(self.tab_Output, orient="horizontal")
        separator.grid(column=0, row=2, columnspan=4, sticky=W + E, padx=self.padx)
        self.tab_Output.columnconfigure(3, weight=1)

        arr_out = [
            "Create Multiple Input Files",
            "Number of Ini Files",
            "Pertubutions-Population(min,max)",
            "Pertubutions-Generation(min,max)",
            "Pertubutions-Mutation(min,max)",
        ]
        self.description_tabs(arr_out, self.tab_Output, row=[3, 5, 6, 7, 8])

        # Create New pertubutuions
        checkbutton_pertub = ttk.Checkbutton(
            self.tab_Output, var=pertub_check, command=checkbox_multi
        )
        checkbutton_pertub.grid(column=1, row=3, sticky=W + E, padx=self.padx)

        pertub_list = list(range(1, 101))

        text = "Each entry allows user to control perturbation percentage of the desire variables."
        entry = ttk.Label(self.tab_Output, text=text, font=self.labelFont)
        entry.grid_configure(
            column=0, row=4, columnspan=3, sticky=W + E, padx=self.padx, pady=self.pady
        )

        entry_n_ini = tk.Entry(
            self.tab_Output, textvariable=self.n_ini, font=self.entryFont
        )
        entry_n_ini.grid(column=1, row=5, columnspan=2, sticky=(W, E), padx=self.padx)

        entry_n_ini.config(state="disabled")

        width = 5
        # --------------
        entry_pertub_pop_min = ttk.Entry(
            self.tab_Output, width=width, textvariable=pop_min, font=self.entryFont
        )
        entry_pertub_pop_min.grid(column=1, row=6, sticky=(W, E), padx=self.padx)

        entry_pertub_pop_max = ttk.Entry(
            self.tab_Output, width=width, textvariable=pop_max, font=self.entryFont
        )
        entry_pertub_pop_max.grid(column=2, row=6, sticky=(W, E), padx=self.padx)

        entry_pertub_pop_min.config(state="disabled")
        entry_pertub_pop_max.config(state="disabled")

        # --------------
        entry_pertub_gen_min = ttk.Entry(
            self.tab_Output, width=width, textvariable=gen_min, font=self.entryFont
        )
        entry_pertub_gen_min.grid(column=1, row=7, sticky=(W, E), padx=self.padx)

        entry_pertub_gen_max = ttk.Entry(
            self.tab_Output, width=width, textvariable=gen_max, font=self.entryFont
        )
        entry_pertub_gen_max.grid(column=2, row=7, sticky=(W, E), padx=self.padx)

        entry_pertub_gen_min.config(state="disabled")
        entry_pertub_gen_max.config(state="disabled")

        # --------------
        entry_pertub_mut_min = ttk.Entry(
            self.tab_Output, width=width, textvariable=mut_min, font=self.entryFont
        )
        entry_pertub_mut_min.grid(column=1, row=8, sticky=(W, E), padx=self.padx)

        entry_pertub_mut_max = ttk.Entry(
            self.tab_Output, width=width, textvariable=mut_max, font=self.entryFont
        )
        entry_pertub_mut_max.grid(column=2, row=8, sticky=(W, E), padx=self.padx)

        entry_pertub_mut_min.config(state="disabled")
        entry_pertub_mut_max.config(state="disabled")

        # --------------

        button_gen_nini = tk.Button(
            self.tab_Output, text="Generate Input Files", command=generate_multi_ini
        )
        button_gen_nini.grid(
            column=0, row=9, columnspan=3, sticky=W + E, padx=self.padx, pady=self.pady
        )
        button_gen_nini.config(state="disabled")

        button_run_nini = tk.Button(
            self.tab_Output, text="Run Multiple Input Files", command=run_multi_ini
        )
        button_run_nini.grid(
            column=0, row=10, columnspan=3, sticky=W + E, padx=self.padx, pady=self.pady
        )
        button_run_nini.config(state="disabled")

    def Build_analysis_tab(self):

        def select_analysis_folder():
            folder_name = filedialog.askdirectory(
                initialdir=os.getcwd(), title="Select folder"
            )
            if not folder_name:
                analysis_folder.set("Please choose a directory")
            else:
                # folder_name = os.path.join(folder_name,'feff')
                analysis_folder.set(folder_name)
                # print(self.feff_file.get())

        def run_analysis():
            params = {}
            # params['base'] = Path(os.getcwd()).parent
            params["base"] = ""
            params["Kmin"] = self.kmin.get()
            params["Kmax"] = self.kmax.get()
            params["kweight"] = self.k_weight.get()
            params["deltak"] = self.delta_k.get()
            params["rbkg"] = self.r_bkg.get()
            params["bkgkw"] = self.bkg_kw.get()
            params["bkgkmax"] = self.bkg_kmax.get()
            params["front"] = self.feff_file.get().split(",")
            params["CSV"] = self.data_file.get()
            params["optimize"] = self.path_optimize.get()
            params["series_index"] = series_index.get()
            params["series"] = self.series.get()
            # set up the params
            self.txtbox.insert(tk.END, "Running Analysis...\n")
            analysis_plot.setup_params(params)

            paths = "[" + self.path_list.get() + "]"
            paths_list = ast.literal_eval(paths)
            # analysis_plot.setup_paths(paths)
            print(paths)
            print(paths_list)
            print(type(paths_list))
            analysis_plot.setup_paths(ast.literal_eval(paths))

            analysis_plot.setup_dirs(analysis_folder.get())
            arr_str, latex_str = analysis_plot.extract_and_run(analysis_folder)

            self.txtbox.insert(tk.END, arr_str)
            self.txtbox.insert(tk.END, "-------------\n")
            self.txtbox.insert(tk.END, latex_str)
            self.txtbox.insert(tk.END, "\n-------------\n")
            self.txtbox.insert(tk.END, "Done")

        def plot_occurances():
            # plot the self_occurances, and then plot the others
            """
            To do: occurances only works with one
            """
            analysis_plot.plot_occurrences(
                analysis_folder.get(),
                self.path_optimize_pert.get(),
                self.path_list.get(),
            )

        analysis_folder = StringVar(self.tab_Analysis, "Please choose a directory")
        series_index = IntVar(self.tab_Analysis, 0)
        self.tab_Analysis.columnconfigure(1, weight=1)
        self.tab_Analysis.columnconfigure(2, weight=1)
        self.tab_Analysis.rowconfigure(10, weight=1)

        arr_out = ["Select folder"]
        self.description_tabs(arr_out, self.tab_Analysis)

        entry_analysis_folder = tk.Entry(
            self.tab_Analysis, textvariable=analysis_folder, font=self.entryFont
        )
        entry_analysis_folder.grid(
            column=1, row=0, columnspan=2, sticky=W + E, padx=self.padx, pady=self.pady
        )

        button_analysis_folder = ttk.Button(
            self.tab_Analysis,
            text="Choose",
            command=select_analysis_folder,
            style="my.TButton",
        )
        button_analysis_folder.grid(
            column=3, row=0, sticky=W + E, padx=self.padx, pady=self.pady
        )

        separator = ttk.Separator(self.tab_Analysis, orient="horizontal")
        separator.grid(column=0, row=1, columnspan=4, sticky=W + E, padx=self.padx)

        analysis_plot = Analysis_Plot(self.tab_Analysis, self.mylarch)

        # ----------------------------------------------------------------------
        button_Run_Analysis = ttk.Button(
            self.tab_Analysis,
            text="Run Analysis",
            command=run_analysis,
            style="my.TButton",
        )
        button_Run_Analysis.grid(
            column=0, row=2, columnspan=4, sticky=W + E, padx=self.padx, pady=self.pady
        )

        button_plot_kR = ttk.Button(
            self.tab_Analysis,
            text="Plot K and R Spectrum",
            command=analysis_plot.plot_k_r_space,
            style="my.TButton",
        )
        button_plot_kR.grid(
            column=0, row=3, columnspan=1, sticky=W + E, padx=self.padx, pady=self.pady
        )

        button_plot_individual = ttk.Button(
            self.tab_Analysis,
            text="Plot K Individual",
            command=analysis_plot.plot_individual,
            style="my.TButton",
        )
        button_plot_individual.grid(
            column=1, row=3, columnspan=1, sticky=W + E, padx=self.padx, pady=self.pady
        )

        button_plot_error = ttk.Button(
            self.tab_Analysis,
            text="Plot Error",
            command=analysis_plot.plot_error,
            style="my.TButton",
        )
        button_plot_error.grid(
            column=2, row=3, columnspan=1, sticky=W + E, padx=self.padx, pady=self.pady
        )

        button_plot_occurances = ttk.Button(
            self.tab_Analysis,
            text="Plot Occurances",
            command=plot_occurances,
            style="my.TButton",
        )
        button_plot_occurances.grid(
            column=3, row=3, columnspan=1, sticky=W + E, padx=self.padx, pady=self.pady
        )

    def Build_expert_tab(self):

        def select_sabcor_file():
            file_name = filedialog.askopenfilename(
                initialdir=os.getcwd(),
                title="Choose Sabcor Input File",
                filetypes=(("inp file", "*.inp"), ("all files", "*.*")),
            )
            if not file_name:
                self.sabcor_input_file.set("Please choose a file")
            else:
                self.sabcor_input_file.set(file_name)

        def select_sabcor_executable():
            file_name = filedialog.askopenfilename(
                initialdir=os.getcwd(),
                title="Choose Sabcor Executable",
                filetypes=(("exe file", "*.exe"), ("all files", "*.*")),
            )
            if not file_name:
                self.sabcor_executable.set("Please choose a file")
            else:
                self.sabcor_executable.set(file_name)

        def checkbox_sabcor():
            if self.sabcor_toggle.get() == True:
                sabcor_input_file.config(state="normal")
                button_sabcor_input_file.config(state="normal")
                sabcor_executable.config(state="normal")
                button_sabcor_executable.config(state="normal")
                # allow only single entry
                self.ncomp.set(1)
                entry_ncomp.config(state="disabled")

            else:
                sabcor_input_file.config(state="disabled")
                button_sabcor_input_file.config(state="disabled")
                sabcor_executable.config(state="disabled")
                button_sabcor_executable.config(state="disabled")
                entry_ncomp.config(state="normal")

        arr_expert = [
            "Selector Options",
            "Crossover Options",
            "Solver Options",
            "Override Num Compounds",
            "Sabcor Toggle",
            "Sabcor Input File",
            "Sabcor Executable",
        ]
        self.description_tabs(arr_expert, self.tab_Expert, row=[0, 1, 2, 3, 5, 6, 7])

        self.tab_Expert.grid_columnconfigure(1, weight=1)

        sel_list = ["Roulette Wheel", "Tournament"]
        selection_combo = ttk.Combobox(
            self.tab_Expert,
            width=7,
            values=sel_list,
            textvariable=self.selection_options,
            font=self.entryFont,
        )
        selection_combo.grid(
            column=1, row=0, columnspan=3, sticky=(W, E), padx=self.padx
        )

        cro_list = [
            "Uniform Crossover",
            "Single Point Corssover",
            "Dual Point Crossover",
            "Arithmetic Crossover",
            "Or Crossover",
            "Average Crossover",
        ]
        crossover_combo = ttk.Combobox(
            self.tab_Expert,
            width=7,
            values=cro_list,
            textvariable=self.crossover_options,
            font=self.entryFont,
        )
        crossover_combo.grid(
            column=1, row=1, columnspan=3, sticky=(W, E), padx=self.padx
        )

        sol_list = [
            "Genetic Algorithm",
            "Genetic Algorithm with Rechenberg",
            "Differential Evolution",
        ]
        solver_combo = ttk.Combobox(
            self.tab_Expert,
            width=7,
            values=sol_list,
            textvariable=self.solver_options,
            font=self.entryFont,
        )
        solver_combo.grid(column=1, row=2, columnspan=3, sticky=(W, E), padx=self.padx)

        ncomp_list = list(range(1, 101))
        entry_ncomp = ttk.Combobox(
            self.tab_Expert,
            width=7,
            values=ncomp_list,
            textvariable=self.ncomp,
            font=self.entryFont,
        )
        entry_ncomp.grid(column=1, row=3, columnspan=3, sticky=(W, E), padx=self.padx)

        separator = ttk.Separator(self.tab_Expert, orient="horizontal")
        separator.grid(column=0, row=4, columnspan=4, sticky=W + E, padx=self.padx)

        sabcor_toggle = ttk.Checkbutton(
            self.tab_Expert, variable=self.sabcor_toggle, command=checkbox_sabcor
        )
        sabcor_toggle.grid(column=1, row=5, sticky=(W, E), padx=self.padx)
        if self.os == "Windows":
            sabcor_toggle.configure(state="disabled")

        sabcor_input_file = tk.Entry(
            self.tab_Expert, textvariable=self.sabcor_input_file, font=self.entryFont
        )
        sabcor_input_file.grid(
            column=1, row=6, sticky=(W, E), padx=self.padx, pady=self.pady
        )
        sabcor_input_file.config(state="disabled")

        button_sabcor_input_file = ttk.Button(
            self.tab_Expert,
            text="Choose",
            command=select_sabcor_file,
            style="my.TButton",
        )
        button_sabcor_input_file.grid(
            column=3, row=6, sticky=W, padx=self.padx, pady=self.pady
        )
        button_sabcor_input_file.config(state="disabled")

        sabcor_executable = tk.Entry(
            self.tab_Expert, textvariable=self.sabcor_executable, font=self.entryFont
        )
        sabcor_executable.grid(
            column=1, row=7, sticky=(W, E), padx=self.padx, pady=self.pady
        )
        sabcor_executable.config(state="disabled")

        button_sabcor_executable = ttk.Button(
            self.tab_Expert,
            text="Choose",
            command=select_sabcor_executable,
            style="my.TButton",
        )
        button_sabcor_executable.grid(
            column=3, row=7, sticky=W, padx=self.padx, pady=self.pady
        )
        button_sabcor_executable.config(state="disabled")

    def On_closing(self):
        """
        on closing function
        """
        if messagebox.askokcancel("Quit", "Do you want to quit?"):
            self.stop_term()
            if hasattr(self, "terminal"):
                self.root.quit()
                self.terminal.destroy()
            else:
                self.root.quit()

    def Run(self):
        """
        Run the code
        """
        self.root.protocol("WM_DELETE_WINDOW", self.On_closing)
        self.beta_popup()

        self.root.mainloop()

    def beta_popup(self):
        beta_popup = tk.Toplevel(self.root)
        beta_popup.wm_title("Warning")
        msg = f"This Graphical User Interface is still under active development.\nPlease contact us using Github Issues.\n EXAFS Neo: {__version__}"
        entry = ttk.Label(beta_popup, text=msg, justify="center")
        entry.grid(column=0, row=0, padx=5, pady=3)
        B1 = ttk.Button(beta_popup, text="Okay", command=beta_popup.destroy)
        B1.grid(column=0, row=1, padx=5, pady=3)

        beta_popup.grid_columnconfigure((0, 1), weight=1)
        beta_popup.grid_rowconfigure((0, 1), weight=1)
        beta_popup.protocol("WM_DELETE_WINDOW", beta_popup.destroy)
        beta_popup.attributes("-topmost", "true")


def main():
    root = App()
    root.Run()


if __name__ == "__main__":
    main()
