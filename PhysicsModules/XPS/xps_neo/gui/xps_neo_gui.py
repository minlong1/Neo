"""
Authors     Alaina Humiston, Megan Burrill, Miu lun lau
Email       athompson9@hawk.iit.edu, mburrill@hawk.iit.edu, andylau@u.boisestate.edu
Version     0.2
Date        Mar 15, 2024
"""

"""
TODO
- analysis [Done]
- improve graphing
- preprocessing?
- connect calibration to nano-indent
"""
#Testing to see if I can push with other computer
from tkinter import *
from threading import Thread
from tkinter import ttk, Tk, N, W, E, S, StringVar, IntVar, DoubleVar, BooleanVar, Checkbutton, NORMAL, DISABLED, \
    scrolledtext, filedialog, messagebox, LabelFrame, Toplevel, END, TOP
from tkinter.font import Font
from tokenize import Double

# import matplotlib
import matplotlib
from matplotlib.hatch import HorizontalHatch
import numpy as np
import configparser

import time


matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
# from threading import *
# import os
import os, subprocess, asyncio
import signal
import pathlib
import multiprocessing as mp
# custom Libraries
from PhysicsModules.XPS.xps_neo.gui.xps_plot import Data_plot, Analysis_plot
#import preprocess_data
from PhysicsModules.XPS.xps_neo import xps_data as data
from PhysicsModules.XPS.xps_neo.gui import xps_analysis2


def _check_cli_version():
    """Warn if the `xps_neo` CLI on PATH is not this package's version.

    The GUI runs fits by spawning the installed CLI; a stale install
    silently fits with old code while the GUI analyzes with new code
    (Phase 5 deployment guard - see CLAUDE.md 'reinstall gotcha').
    """
    import xps_neo
    try:
        out = subprocess.run(["xps_neo", "--version"], capture_output=True,
                             text=True, timeout=15).stdout.strip()
    except Exception as exc:
        print(f"WARNING: could not run the xps_neo CLI ({exc}); "
              f"fits will not work until it is installed (pip install -e .)")
        return
    expected = f"xps_neo {xps_neo.__version__}"
    if out != expected:
        print(f"WARNING: the xps_neo CLI on PATH reports {out!r} but the "
              f"GUI imports package version {xps_neo.__version__}. Fits and "
              f"analysis would use different code - reinstall with "
              f"'pip install -e .' in the matching environment.")


_check_cli_version()
from PhysicsModules.XPS.xps_neo.periodic_table import ElementData
#from uncertainties import ufloat


#-----------------TO DO LIST-----------------
#
# Fix Stop button in code --> This may be a Windows issue. Look into proc.kill
#
# Shirley (Not Shirley-Sherwood) background has issue with mutiple backgrounds selected: Hard coded to make work in xps.py --> Need to actually fix
#
# Error if background parameter range is too small: Test and fix
#
# Shirley Background: Need to unrestrain on low BE side of data to allow for large Lorentzian values --> What to use as reference? Baseline? Always have baseline selected and use value inside Shirley formula?
#
# Change name of parameters for less confusion: GFWHM and LFWHM. Too much mixing of Gaussian with sigma and Lorentzian with gamma
#
# Find a way to disable Doublet buttons. Also Coster-Kronig if gaussian or lorentzian or doniach-sunjic peak shape is selected
#
# Check for if fit is low on low BE side (around BE + LFWHM): Increase LFWHM --> If issue is still there then consider adding a peak. How many gens? 25, 50?
#
#Fix convolve method: Makes points on high BE side of data tip down. Issue in scipy convolution. Manual convolution?
#
# Python 3.14 not compatible? 3.12 works though
#
# Change .fil for aanalyzer to incorporate Coster-Kronig peaks








class App():
    """
    Start of the application
    """

    def __init__(self):
        self.__version__ = 0.1
        self.root = Tk(className='XPS Neo GUI')
        self.root.wm_title("XPS GUI (Pre-Beta)")
        self.root.geometry("1000x650") #Changed size from 975x650 to 1100x750
        self.mainframe = ttk.Notebook(self.root, padding='3')
        self.mainframe.grid(column=0, row=0, sticky=(N, W, E, S), columnspan=5)
        self.root.columnconfigure(0, weight=1)
        self.root.columnconfigure(1, weight=1)
        self.root.rowconfigure(0, weight=1)
        self.root.resizable(True, True)
        self.padx = 2
        self.pady = 1

        # Specify standard font
        self.entryFont = Font(family='TkTextFont', size=10) #Changed from 11
        self.labelFont = Font(family='TkMenuFont', size=10) # Changed from 11
        s = ttk.Style()
        # s.configure('my.TButton', font=labelFont)

        # Set multiprocessing run type
        mp.set_start_method('spawn')

        # initialize variables
        self.initialize_var()
        self.initialize_tabs()
        self.build_tabs()


    



    def initialize_var(self):
        """
        Initalize all possible variables in the GUI
        """
        # Inputs (column 0)
        # Averaging Section of Input Tab:
        self.all_files_display = StringVar(self.root, 'Please choose file/s')  #Files that you want to average
        self.all_files = np.array([])
        self.averaged_file = StringVar(self.root, 'Please choose file/s') #File to save the averaged data to
        self.num_files = DoubleVar(self.root, 0) #Number of files to be averaged after selection
        self.num_points = DoubleVar(self.root, 0) #Number of points in each data file
        self.xPoints_for_data_plot = np.array([])
        self.yPoints_for_data_plot = np.array([])
        self.xPoints_avg = np.array([]) #Averaged x points for plotting
        self.yPoints_avg = np.array([]) #Averaged y points for plotting
        self.plot_data = np.array([],[])

        # There is a string var to give the prompt on selection, then a path saved to allow easier and safer manipulation
        self.csv_file = StringVar(self.root, "Please choose a data file")
        self.csv_generate_from = pathlib.Path()
        self.file_name = "Data file"
        self.csv_folder = StringVar(self.root, "Please choose a folder containing data files")
        self.csv_folder_path = pathlib.Path()
        self.output_folder = StringVar(self.root, 'Please choose a folder to save outputs')
        self.output_folder_path = pathlib.Path()
        self.output_file = pathlib.Path()
        #File reading options
        self.skipLn = StringVar(self.root,"0")
        self.x_offset = StringVar(self.root,"0")
        self.y_offset = StringVar(self.root,"0")
        self.scale_var = BooleanVar(self.root,"False")

        

        #Adding in fit file selection to input variables used in previous fitting process
        self.fit_file = StringVar(self.root, "Load in fit file information")
        self.fit_generate_from = pathlib.Path()
        self.fit_file_name = "Data file"
        self.fit_file_selected = False
        #self.fit_peak_number = 1

        self.poly1_value = ''
        self.poly2_value = ''
        self.poly3_value = ''
        self.toug2_value = ''
        self.toug3_value = ''
        self.linear_value = ''
        self.exp_value = ''
        self.shirley_value = ''
        self.reg_shirley_value = ''
        


        # Waiting to create path variable for calibration because I think it will be deleted
        #self.csv_calibration_file = StringVar(self.root, "Please choose a calibration file") #--------------Will need to comment out-----------
        self.yes_folder = IntVar()
        self.multi_known = BooleanVar(value =False)
        self.multi_known = False  # This will change to true when the user has already said if they do/do not want to generate/run multiple files in the folder so generate ini does nto keep asking
        self.filelist = []
        self.command_list = mp.Queue()
        self.proc_list = []
        self.pid_list = []
        # Preprocessing
        self.preprocess_file = pathlib.Path()
        self.stop_not_pressed = True

        # Variables for the dropdown menu
        self.file_menu = ttk.Combobox()
        self.data_obj = 0

        self.noise_check = True
        self.range_check = True

        # Populations (column 1)
        self.population = IntVar(self.root, 2000)
        self.num_gen = IntVar(self.root, 10)
        self.best_sample = IntVar(self.root, 20)
        self.lucky_few = IntVar(self.root, 10)

        # Mutations (column 2)
        self.chance_of_mutation = IntVar(self.root, 20)
        self.original_chance_of_mutation = IntVar(self.root, 20) #I dont think this is being used anywhere
        self.mutated_options = IntVar(self.root, 0)

        # Fitting parameters(column 4)

        #New XPS parameters
        #how to make unique for each peak?
        self.number_of_peaks = DoubleVar(self.root, 1)
        self.element_select = []
        self.photoLine_select = []
        self.BE_guesses = []
        self.sigma_guesses = []
        self.gamma_guesses = []
        self.amp_guesses = []
        self.asym_guesses = []
        self.sos_guesses = []
        self.br_guesses = []
        #self.so_split = []

        # Parameter ranges (column 5)
        self.BE_up_lim = []
        self.BE_low_lim = []
        self.BE_limit = []
        self.BE_corr = []
        self.BE_corr_mult = []

        self.sigma_up_lim = []
        self.sigma_low_lim = []
        self.sigma_limit = []
        self.sigma_corr = []
        self.sigma_corr_mult = []

        self.gamma_up_lim = []
        self.gamma_low_lim = []
        self.gamma_limit = []
        self.gamma_corr = []
        self.gamma_corr_mult = []

        self.amp_up_lim = []
        self.amp_low_lim = []
        self.amp_limit = []
        self.amp_corr = []
        self.amp_corr_mult = []

        self.asym_up_lim = []
        self.asym_low_lim = []
        self.asym_limit = []
        self.asym_corr = []
        self.asym_corr_mult = []

        self.sos_up_lim = []
        self.sos_low_lim = []
        self.sos_limit = []
        self.sos_corr = []
        self.sos_corr_mult = []

        self.br_up_lim = []
        self.br_low_lim = []
        self.br_limit = []
        self.br_corr = []
        self.br_corr_mult = []

       


        for i in range(10):
            self.sigma_guesses.append(DoubleVar(self.root, 1.0))
            self.gamma_guesses.append(DoubleVar(self.root, 0.25))
            self.amp_guesses.append(DoubleVar(self.root, 2000))
            self.BE_guesses.append(DoubleVar(self.root, 284.50))
            #self.so_split.append(DoubleVar(self.root, 0.00))
            self.sos_guesses.append(DoubleVar(self.root, 0.00))
            self.br_guesses.append(DoubleVar(self.root, 0.5))
            self.asym_guesses.append(DoubleVar(self.root, 1)) #Assume zero asymmetry at first
            
            self.BE_up_lim.append(DoubleVar(self.root, 0.2))
            self.BE_low_lim.append(DoubleVar(self.root, -0.2))
            self.sigma_up_lim.append(DoubleVar(self.root, 0.25))
            self.sigma_low_lim.append(DoubleVar(self.root, -0.25))
            self.gamma_up_lim.append(DoubleVar(self.root, 0.13))
            self.gamma_low_lim.append(DoubleVar(self.root, -0.13))
            self.amp_up_lim.append(DoubleVar(self.root, 500))
            self.amp_low_lim.append(DoubleVar(self.root, -500))
            self.asym_up_lim.append(DoubleVar(self.root, 10))
            self.asym_low_lim.append(DoubleVar(self.root, 1))
            self.sos_up_lim.append(DoubleVar(self.root, 0.2))
            self.sos_low_lim.append(DoubleVar(self.root, -0.2))
            self.br_up_lim.append(DoubleVar(self.root, 0.01))
            self.br_low_lim.append(DoubleVar(self.root, -0.01))
           

        self.peaks = []
        self.path_branching = []
        self.peak_singlet = []
        self.peak_coster_kronig = []
        self.data_KE = False
        
        self.data_XES = False
        self.data_PT = False
        self.LFWHM_alt = BooleanVar(self.root, False)
        self.data_peak_add = BooleanVar(self.root, False)
        for i in range(10):
            self.peaks.append(StringVar(self.root, "Select a peak type"))
            self.element_select.append(StringVar(self.root, " "))
            self.photoLine_select.append(StringVar(self.root, " "))
            self.path_branching.append(DoubleVar(self.root, 0.5))
            self.peak_singlet.append(BooleanVar(self.root, True))
            self.peak_coster_kronig.append(BooleanVar(self.root, False))
            


            self.BE_limit.append(BooleanVar(self.root, False))
            self.BE_corr.append(StringVar(self.root, "Peak #"))
            self.BE_corr_mult.append(DoubleVar(self.root, 1.0))

            self.sigma_limit.append(BooleanVar(self.root, False))
            self.sigma_corr.append(StringVar(self.root, "Peak #"))
            self.sigma_corr_mult.append(DoubleVar(self.root, 1.0))

            self.gamma_limit.append(BooleanVar(self.root, False))
            self.gamma_corr.append(StringVar(self.root, "Peak #"))
            self.gamma_corr_mult.append(DoubleVar(self.root, 1.0))

            self.amp_limit.append(BooleanVar(self.root, False))
            self.amp_corr.append(StringVar(self.root, "Peak #"))
            self.amp_corr_mult.append(DoubleVar(self.root, 1.0))

            self.asym_limit.append(BooleanVar(self.root, False))
            self.asym_corr.append(StringVar(self.root, "Peak #"))
            self.asym_corr_mult.append(DoubleVar(self.root, 1.0))

            self.sos_limit.append(BooleanVar(self.root, False))
            self.sos_corr.append(StringVar(self.root, "Peak #"))
            self.sos_corr_mult.append(DoubleVar(self.root, 1.0))

            self.br_limit.append(BooleanVar(self.root, False))
            self.br_corr.append(StringVar(self.root, "Peak #"))
            self.br_corr_mult.append(DoubleVar(self.root, 1.0))


            #self.data_KE.append(BooleanVar(self.root, False))
        self.background_types = []
        #self.path_bkgn = StringVar(self.root, "Select a background type") #Dont think we need anymore because we want to use an array now
        #self.path_branching = StringVar(self.root, "Select branching ratio")

        #self.spinOrbitSplit_min, self.spinOrbitSplit_max, self.spinOrbitSplit_delta = -0.2, 0.2, 0.01
        #self.branching_ratio_min, self.branching_ratio_max, self.branching_ratio_delta = -0.02, 0.02, 0.01
        self.gamma_CK_min, self.gamma_CK_max, self.gamma_CK_delta = -0.01, 2.0, 0.0001 #What should the upper limit be? x2 lorentz value?
        self.asym_CK_min, self.asym_CK_max, self.asym_CK_delta = 1.0, 20.0, 0.001
        self.BE_min,self.BE_max,self.BE_delta = -0.5,0.5,.01
        #self.sigma_min,self.sigma_max,self.sigma_delta = 0,2,.001
        #self.fwhm_min,self.fwhm_max,self.fwhm_delta = 0,0.5,.001
        #self.amp_min,self.amp_max,self.amp_delta = 0,5000, 0.05 #Was 0-5000: arbitrary. Max set in xps.py

        #GOING TO GET RID OF THESE PARAMETERS: INSTEAD HAVE THE PARAMETER RANGES TAB
        
        self.sigma_min,self.sigma_max,self.sigma_delta = -0.13,0.13,0.0001 #changed from +/- 0.13
        self.fwhm_min,self.fwhm_max,self.fwhm_delta = -0.13,0.13,0.0001
        self.amp_min,self.amp_max,self.amp_delta = -500,500, 0.05
        self.asymmetry_min,self.asymmetry_max,self.asymmetry_delta = 1.0,15.0, 0.001 #NEED TO ALLOW USER TO INPUT GUESS
        self.asymmetryDoniach_min,self.asymmetryDoniach_max,self.asymmetryDoniach_delta = 0, 1, 0.001 #Doniach Sunjic asymmetry
        self.background_min,self.background_max,self.background_delta = 0,5000,0.05 #Max set in xps.py
        self.slope_min,self.slope_max,self.slope_delta = 0,200, 0.0001 #Max AND Min set in xps.py
        self.slope_min,self.slope_max,self.slope_delta = 0,200, 0.0001 #Max AND Min set in xps.py
        

        # Graph (column 6)
        self.print_graph = BooleanVar(self.root, False)
        self.num_output_paths = BooleanVar(self.root, True)
        self.steady_state_exit = BooleanVar(self.root, True)

        # Expert (column 7)
        self.expertLn = StringVar(self.root,"0")
        self.expert = False
        self.gen_late_stage = IntVar(self.root, 1000)


        # Output tab (column 8)
        self.print_graph = BooleanVar(self.root, False)
        self.num_output_paths = BooleanVar(self.root, True)
        self.steady_state_exit = BooleanVar(self.root, True)
        self.n_ini = IntVar(self.root, 100)
        self.pop_min = IntVar(self.root, 100)
        self.pop_max = IntVar(self.root, 5001)
        self.gen_min = IntVar(self.root, 20)
        self.gen_max = IntVar(self.root, 501)
        self.mut_min = IntVar(self.root, 20)
        self.mut_max = IntVar(self.root, 51)
        self.run_folder = BooleanVar(self.root, False)
        self.pertub_check = IntVar(self.root, 0)
        self.checkbutton_whole_folder = ttk.Checkbutton()

        self.analysis_dir = StringVar(self.root, "Please choose a data file")
        # Analysis (column 7)
        # I don't know what goes here so leave blank instead


        #Added this code in to update the amplitude ranges in the param_range_tab so that when you update the amplitude value for any of the peaks, it resets that amplitude upper and lower range to be within reasonable values. 
        self.number_of_peaks


        def amp_change_peak_1(*args):
            peak_var = 0
            try: 
                   
                self.amp_up_lim[peak_var] = DoubleVar(self.root, self.amp_guesses[peak_var].get()*0.2) 
                self.amp_low_lim[peak_var] = DoubleVar(self.root, -self.amp_guesses[peak_var].get()) 
            except:
                self.amp_up_lim[peak_var] = 0.0
                self.amp_low_lim[peak_var] = 0.0
                                        
            self.build_param_range_tab()
        
        def amp_change_peak_2(*args):
            peak_var = 1
            try: 
                   
                self.amp_up_lim[peak_var] = DoubleVar(self.root, self.amp_guesses[peak_var].get()*0.2) 
                self.amp_low_lim[peak_var] = DoubleVar(self.root, -self.amp_guesses[peak_var].get()) 
            except:
                self.amp_up_lim[peak_var] = 0.0
                self.amp_low_lim[peak_var] = 0.0
                                        
            self.build_param_range_tab()
        
        def amp_change_peak_3(*args):
            peak_var = 2
            try: 
                   
                self.amp_up_lim[peak_var] = DoubleVar(self.root, self.amp_guesses[peak_var].get()*0.2) 
                self.amp_low_lim[peak_var] = DoubleVar(self.root, -self.amp_guesses[peak_var].get()) 
            except:
                self.amp_up_lim[peak_var] = 0.0
                self.amp_low_lim[peak_var] = 0.0
                                        
            self.build_param_range_tab()

        def amp_change_peak_4(*args):
            peak_var = 3
            try: 
                   
                self.amp_up_lim[peak_var] = DoubleVar(self.root, self.amp_guesses[peak_var].get()*0.2) 
                self.amp_low_lim[peak_var] = DoubleVar(self.root, -self.amp_guesses[peak_var].get()) 
            except:
                self.amp_up_lim[peak_var] = 0.0
                self.amp_low_lim[peak_var] = 0.0
                                        
            self.build_param_range_tab()

        def amp_change_peak_5(*args):
            peak_var = 4
            try: 
                   
                self.amp_up_lim[peak_var] = DoubleVar(self.root, self.amp_guesses[peak_var].get()*0.2) 
                self.amp_low_lim[peak_var] = DoubleVar(self.root, -self.amp_guesses[peak_var].get()) 
            except:
                self.amp_up_lim[peak_var] = 0.0
                self.amp_low_lim[peak_var] = 0.0
                                        
            self.build_param_range_tab()

        def amp_change_peak_6(*args):
            peak_var = 5
            try: 
                   
                self.amp_up_lim[peak_var] = DoubleVar(self.root, self.amp_guesses[peak_var].get()*0.2) 
                self.amp_low_lim[peak_var] = DoubleVar(self.root, -self.amp_guesses[peak_var].get()) 
            except:
                self.amp_up_lim[peak_var] = 0.0
                self.amp_low_lim[peak_var] = 0.0
                                        
            self.build_param_range_tab()

        def amp_change_peak_7(*args):
            peak_var = 6
            try: 
                   
                self.amp_up_lim[peak_var] = DoubleVar(self.root, self.amp_guesses[peak_var].get()*0.2) 
                self.amp_low_lim[peak_var] = DoubleVar(self.root, -self.amp_guesses[peak_var].get()) 
            except:
                self.amp_up_lim[peak_var] = 0.0
                self.amp_low_lim[peak_var] = 0.0
                                        
            self.build_param_range_tab()

        def amp_change_peak_8(*args):
            peak_var = 7
            try: 
                   
                self.amp_up_lim[peak_var] = DoubleVar(self.root, self.amp_guesses[peak_var].get()*0.2) 
                self.amp_low_lim[peak_var] = DoubleVar(self.root, -self.amp_guesses[peak_var].get()) 
            except:
                self.amp_up_lim[peak_var] = 0.0
                self.amp_low_lim[peak_var] = 0.0
                                        
            self.build_param_range_tab()

        def amp_change_peak_9(*args):
            peak_var = 8
            try: 
                   
                self.amp_up_lim[peak_var] = DoubleVar(self.root, self.amp_guesses[peak_var].get()*0.2) 
                self.amp_low_lim[peak_var] = DoubleVar(self.root, -self.amp_guesses[peak_var].get()) 
            except:
                self.amp_up_lim[peak_var] = 0.0
                self.amp_low_lim[peak_var] = 0.0
                                        
            self.build_param_range_tab()

        def amp_change_peak_10(*args):
            peak_var = 9
            try: 
                   
                self.amp_up_lim[peak_var] = DoubleVar(self.root, self.amp_guesses[peak_var].get()*0.2) 
                self.amp_low_lim[peak_var] = DoubleVar(self.root, -self.amp_guesses[peak_var].get()) 
            except:
                self.amp_up_lim[peak_var] = 0.0
                self.amp_low_lim[peak_var] = 0.0
                                        
            self.build_param_range_tab()

        amps = {
                0: amp_change_peak_1,
                1: amp_change_peak_2,
                2: amp_change_peak_3,
                3: amp_change_peak_4,
                4: amp_change_peak_5,
                5: amp_change_peak_6,
                6: amp_change_peak_7,
                7: amp_change_peak_8,
                8: amp_change_peak_9,
                9: amp_change_peak_10
            }
        


        def BE_change_peak_1(*args):
            peak_BE_var = 0
            

            try:
                x = []
                y = []
                x.append(self.data_obj.get_x(self.data_KE, self.data_XES))
                y.append(self.data_obj.get_y(self.scale_var))

                close_val = x[0][0]  
                    
                min_diff = abs(x[0][0]- float(self.BE_guesses[peak_BE_var].get())) 
                for k in range(len(x[0])):
                    
                    diff = abs(x[0][k] - float(self.BE_guesses[peak_BE_var].get()))
                    
                    if diff < min_diff:
                        min_diff = diff
                        close_val = k

                close_amp = y[0][close_val]      

                self.amp_guesses[peak_BE_var] = DoubleVar(self.root, close_amp) 
                amps[peak_BE_var]() 
                self.build_fitting_param_tab() #This adds in an extra Baseline value...Gotta get rid of. Also slows down loading significantly when algorithm is running
                self.background_types.remove('Baseline')
                #print("Backgrounds", self.background_types)
            except: 
                pass

        def BE_change_peak_2(*args):
            peak_BE_var = 1

            try:
                x = []
                y = []
                x.append(self.data_obj.get_x(self.data_KE, self.data_XES))
                y.append(self.data_obj.get_y(self.scale_var))

                close_val = x[0][0]  
                    
                min_diff = abs(x[0][0]- float(self.BE_guesses[peak_BE_var].get())) 
                for k in range(len(x[0])):
                    
                    diff = abs(x[0][k] - float(self.BE_guesses[peak_BE_var].get()))
                    
                    if diff < min_diff:
                        min_diff = diff
                        close_val = k

                close_amp = y[0][close_val]      

                self.amp_guesses[peak_BE_var] = DoubleVar(self.root, close_amp) 
                amps[peak_BE_var]() 
                self.build_fitting_param_tab()
                self.background_types.remove('Baseline')
                #print("Backgrounds", self.background_types)

            except: 
                pass

        def BE_change_peak_3(*args):
            peak_BE_var = 2

            try:
                x = []
                y = []
                x.append(self.data_obj.get_x(self.data_KE, self.data_XES))
                y.append(self.data_obj.get_y(self.scale_var))

                close_val = x[0][0]  
                    
                min_diff = abs(x[0][0]- float(self.BE_guesses[peak_BE_var].get())) 
                for k in range(len(x[0])):
                    
                    diff = abs(x[0][k] - float(self.BE_guesses[peak_BE_var].get()))
                    
                    if diff < min_diff:
                        min_diff = diff
                        close_val = k

                close_amp = y[0][close_val]      

                self.amp_guesses[peak_BE_var] = DoubleVar(self.root, close_amp) 
                amps[peak_BE_var]() 
                self.build_fitting_param_tab()
                self.background_types.remove('Baseline')
                #print("Backgrounds", self.background_types)

            except: 
                pass

        def BE_change_peak_4(*args):
            peak_BE_var = 3

            try:
                x = []
                y = []
                x.append(self.data_obj.get_x(self.data_KE, self.data_XES))
                y.append(self.data_obj.get_y(self.scale_var))

                close_val = x[0][0]  
                    
                min_diff = abs(x[0][0]- float(self.BE_guesses[peak_BE_var].get())) 
                for k in range(len(x[0])):
                    
                    diff = abs(x[0][k] - float(self.BE_guesses[peak_BE_var].get()))
                    
                    if diff < min_diff:
                        min_diff = diff
                        close_val = k

                close_amp = y[0][close_val]      

                self.amp_guesses[peak_BE_var] = DoubleVar(self.root, close_amp) 
                amps[peak_BE_var]() 
                self.build_fitting_param_tab()
                self.background_types.remove('Baseline')
                #print("Backgrounds", self.background_types)
            except: 
                pass

        def BE_change_peak_5(*args):
            peak_BE_var = 4

            try:
                x = []
                y = []
                x.append(self.data_obj.get_x(self.data_KE, self.data_XES))
                y.append(self.data_obj.get_y(self.scale_var))

                close_val = x[0][0]  
                    
                min_diff = abs(x[0][0]- float(self.BE_guesses[peak_BE_var].get())) 
                for k in range(len(x[0])):
                    
                    diff = abs(x[0][k] - float(self.BE_guesses[peak_BE_var].get()))
                    
                    if diff < min_diff:
                        min_diff = diff
                        close_val = k

                close_amp = y[0][close_val]      

                self.amp_guesses[peak_BE_var] = DoubleVar(self.root, close_amp) 
                amps[peak_BE_var]() 
                self.build_fitting_param_tab()
                self.background_types.remove('Baseline')
                #print("Backgrounds", self.background_types)

            except: 
                pass

        def BE_change_peak_6(*args):
            peak_BE_var = 5

            try:
                x = []
                y = []
                x.append(self.data_obj.get_x(self.data_KE, self.data_XES))
                y.append(self.data_obj.get_y(self.scale_var))

                close_val = x[0][0]  
                    
                min_diff = abs(x[0][0]- float(self.BE_guesses[peak_BE_var].get())) 
                for k in range(len(x[0])):
                    
                    diff = abs(x[0][k] - float(self.BE_guesses[peak_BE_var].get()))
                    
                    if diff < min_diff:
                        min_diff = diff
                        close_val = k

                close_amp = y[0][close_val]      

                self.amp_guesses[peak_BE_var] = DoubleVar(self.root, close_amp) 
                amps[peak_BE_var]() 
                self.build_fitting_param_tab()
                self.background_types.remove('Baseline')
                #print("Backgrounds", self.background_types)
            except: 
                pass

        def BE_change_peak_7(*args):
            peak_BE_var = 6

            try:
                x = []
                y = []
                x.append(self.data_obj.get_x(self.data_KE, self.data_XES))
                y.append(self.data_obj.get_y(self.scale_var))

                close_val = x[0][0]  
                    
                min_diff = abs(x[0][0]- float(self.BE_guesses[peak_BE_var].get())) 
                for k in range(len(x[0])):
                    
                    diff = abs(x[0][k] - float(self.BE_guesses[peak_BE_var].get()))
                    
                    if diff < min_diff:
                        min_diff = diff
                        close_val = k

                close_amp = y[0][close_val]      

                self.amp_guesses[peak_BE_var] = DoubleVar(self.root, close_amp) 
                amps[peak_BE_var]() 
                self.build_fitting_param_tab()
                self.background_types.remove('Baseline')
                #print("Backgrounds", self.background_types)
            except: 
                pass

        def BE_change_peak_8(*args):
            peak_BE_var = 7

            try:
                x = []
                y = []
                x.append(self.data_obj.get_x(self.data_KE, self.data_XES))
                y.append(self.data_obj.get_y(self.scale_var))

                close_val = x[0][0]  
                    
                min_diff = abs(x[0][0]- float(self.BE_guesses[peak_BE_var].get())) 
                for k in range(len(x[0])):
                    
                    diff = abs(x[0][k] - float(self.BE_guesses[peak_BE_var].get()))
                    
                    if diff < min_diff:
                        min_diff = diff
                        close_val = k

                close_amp = y[0][close_val]      

                self.amp_guesses[peak_BE_var] = DoubleVar(self.root, close_amp) 
                amps[peak_BE_var]() 
                self.build_fitting_param_tab()
                self.background_types.remove('Baseline')
                #print("Backgrounds", self.background_types)
            except: 
                pass

        def BE_change_peak_9(*args):
            peak_BE_var = 8

            try:
                x = []
                y = []
                x.append(self.data_obj.get_x(self.data_KE, self.data_XES))
                y.append(self.data_obj.get_y(self.scale_var))

                close_val = x[0][0]  
                    
                min_diff = abs(x[0][0]- float(self.BE_guesses[peak_BE_var].get())) 
                for k in range(len(x[0])):
                    
                    diff = abs(x[0][k] - float(self.BE_guesses[peak_BE_var].get()))
                    
                    if diff < min_diff:
                        min_diff = diff
                        close_val = k

                close_amp = y[0][close_val]      

                self.amp_guesses[peak_BE_var] = DoubleVar(self.root, close_amp) 
                amps[peak_BE_var]() 
                self.build_fitting_param_tab()
                self.background_types.remove('Baseline')
                #print("Backgrounds", self.background_types)
            except: 
                pass

        def BE_change_peak_10(*args):
            peak_BE_var = 9

            try:
                x = []
                y = []
                x.append(self.data_obj.get_x(self.data_KE, self.data_XES))
                y.append(self.data_obj.get_y(self.scale_var))

                close_val = x[0][0]  
                    
                min_diff = abs(x[0][0]- float(self.BE_guesses[peak_BE_var].get())) 
                for k in range(len(x[0])):
                    
                    diff = abs(x[0][k] - float(self.BE_guesses[peak_BE_var].get()))
                    
                    if diff < min_diff:
                        min_diff = diff
                        close_val = k

                close_amp = y[0][close_val]      

                self.amp_guesses[peak_BE_var] = DoubleVar(self.root, close_amp) 
                amps[peak_BE_var]() 
                self.build_fitting_param_tab()
                self.background_types.remove('Baseline')
                #print("Backgrounds", self.background_types)

            except: 
                pass



        
        #only works for initial input. Took away updating BE from peroidic to make it sitll work with that. Need to figure out why loading in/periodic table updates make this not work
        self.BE_guesses[0].trace_add("write", BE_change_peak_1)
        self.BE_guesses[1].trace_add("write", BE_change_peak_2)
        self.BE_guesses[2].trace_add("write", BE_change_peak_3)
        self.BE_guesses[3].trace_add("write", BE_change_peak_4)
        self.BE_guesses[4].trace_add("write", BE_change_peak_5)
        self.BE_guesses[5].trace_add("write", BE_change_peak_6)
        self.BE_guesses[6].trace_add("write", BE_change_peak_7)
        self.BE_guesses[7].trace_add("write", BE_change_peak_8)
        self.BE_guesses[8].trace_add("write", BE_change_peak_9)
        self.BE_guesses[9].trace_add("write", BE_change_peak_10)
        
        #When reloading parameters this does not update. How to fix?
        self.amp_guesses[0].trace_add("write", amp_change_peak_1)
        self.amp_guesses[1].trace_add("write", amp_change_peak_2)
        self.amp_guesses[2].trace_add("write", amp_change_peak_3)
        self.amp_guesses[3].trace_add("write", amp_change_peak_4)
        self.amp_guesses[4].trace_add("write", amp_change_peak_5)
        self.amp_guesses[5].trace_add("write", amp_change_peak_6)
        self.amp_guesses[6].trace_add("write", amp_change_peak_7)
        self.amp_guesses[7].trace_add("write", amp_change_peak_8)
        self.amp_guesses[8].trace_add("write", amp_change_peak_9)
        self.amp_guesses[9].trace_add("write", amp_change_peak_1)



    def initialize_tabs(self):
        """
        Initialize tabs for the main frame
        """
        s = ttk.Style()
        s.configure('TNotebook.Tab', font=('TkHeadingFont', '11'))
        height = 1
        # Creating tabs
        self.input_tab = ttk.Frame(self.mainframe, height=height)
        self.population_tab = ttk.Frame(self.mainframe, height=height)
        #self.calibration_tab = ttk.Frame(self.mainframe, height=height)
        self.mutation_tab = ttk.Frame(self.mainframe, height=height)
        self.periodicTable_tab = ttk.Frame(self.mainframe, height=height) #New tab. Going to replace Fitting Parameters tab in the future
        self.fitting_param_tab = ttk.Frame(self.mainframe, height=height)
        self.param_range_tab = ttk.Frame(self.mainframe, height=height)
        self.graph_tab = ttk.Frame(self.mainframe, height=height)
        self.expert_tab = ttk.Frame(self.mainframe, height=height)
        self.output_tab = ttk.Frame(self.mainframe, height=height) #Do wo need this tab?
        self.analysis_tab = ttk.Frame(self.mainframe, height=height)

        # Adding tabs
        self.mainframe.add(self.input_tab, text="Inputs")
        self.mainframe.add(self.population_tab, text='Populations')
        self.mainframe.add(self.mutation_tab, text="Mutations")
        self.mainframe.add(self.periodicTable_tab, text="Element Selection") #New tab. Going to replace Fitting Parameters tab in the future
        self.mainframe.add(self.fitting_param_tab, text="Fitting Parameters")
        self.mainframe.add(self.param_range_tab, text= "Parameter Ranges")
        self.mainframe.add(self.graph_tab, text="Plots")
        self.mainframe.add(self.expert_tab, text="Expert")
        self.mainframe.add(self.output_tab, text="Output")
        self.mainframe.add(self.analysis_tab, text="Analysis")

    def build_tabs(self):
        """
        Build tabs. Will call function for each tab
        """
        self.build_global()
        self.build_inputs_tab()
        self.build_population_tab()
        self.build_mutations_tab()
        self.build_fitting_param_tab()
        self.build_param_range_tab()
        self.build_periodicTable_tab()
        self.build_plot_tab()
        self.build_expert_tab()
        self.build_output_tab()
        self.build_analysis_tab()

        self.mainframe.grid_rowconfigure(0, weight=1)
        self.mainframe.grid_columnconfigure(0, weight=1)
        self.mainframe.grid_columnconfigure(1, weight=1)
        self.mainframe.grid_rowconfigure(2, weight=1)
        self.mainframe.grid_columnconfigure(3, weight=1)
        self.mainframe.grid_columnconfigure(4, weight=1)

    def description_tabs(self, arr, tabs, sticky=(W, E), row=None, column=None, return_description=False):
        # Rows = index of rows
        # Loops through array of descriptors to be added to the tabs
        description_list = []
        if row is not None:
            assert len(row) == len(arr)
        for i, inputs in enumerate(arr):
            entry = ttk.Label(tabs, text=inputs, font=self.labelFont)
            if row is not None:
                k = row[i]
            else:
                k = i
            entry.grid_configure(column=column, row=k, sticky=sticky, padx=self.padx, pady=self.pady)
            description_list.append(entry)
        if description_list:
            return description_list

    def description_tabs_column(self, arr, tabs, sticky=(W, E), row=None, column=None, return_description=False):
        # Rows = index of rows
        # Loops through array of descriptors to be added to the tabs
        description_list = []
        if column is not None:
            assert len(column) == len(arr)
        for i, inputs in enumerate(arr):
            entry = ttk.Label(tabs, text=inputs, font=self.labelFont)
            if column is not None:
                k = column[i]
            else:
                k = i
            entry.grid_configure(column=k, row=row, sticky=sticky, padx=self.padx, pady=self.pady)
            description_list.append(entry)
        if description_list:
            return description_list




    # When the user selects a particular file from the input directory in the dropdown menu it will be assigned to the
    # csv generate from variable so that it is used for running
    def file_selected(self, event):
        self.csv_generate_from = pathlib.Path(self.csv_folder_path.joinpath(self.file_menu.get()))
        # print("file selected: ", self.csv_generate_from)

    # Writes the ini file to filename using the user inputs or defaults if nothing is changed
    def write_ini(self, filename):
        # First select data range
        # preprocess_data.read_files(self.csv_generate_from, limits=(self.percent_min.get(), self.percent_max.get()))
        # self.preprocess_file = pathlib.Path.cwd().joinpath('example.txt')

        # print("write ini csv path generate from", str(self.csv_generate_from))
        inputs = ("[Inputs]\ndata_file = {data}\noutput_file = {out} \nskipln = {skipLn} \nx_offset = {x_offset} \ny_offset = {y_offset}"
                  "".format(
            data=str(self.csv_generate_from),
            out=str(self.output_file),
            skipLn = str(self.skipLn.get()),
            x_offset = str(self.x_offset.get()),
            y_offset = str(self.y_offset.get()),
            
        
            #calibration=str(self.csv_calibration_file.get()),#deleted: calibration=str(self.csv_calibration_file.get())
            #dat_cutoff=", ".join(str(i) for i in [self.percent_min.get(), self.percent_max.get()])
        ))
        populations = ("\n\n[Populations]\npopulation = {pop}\nnum_gen = {numgen}\nbest_sample = {best} "
                       "\nlucky_few = {luck}".format(pop=str(self.population.get()),
                                                     numgen=str(self.num_gen.get()),
                                                     best=str(self.best_sample.get()),
                                                     luck=str(self.lucky_few.get())))
        # Sends range for power law equation as [min, max]

        elements = []
        photoelectronLines = []
        guesses = []
        BE_guess_min = []
        BE_guess_max = []
        BE_guess_delta = []
        BE_guesses_range = []
        sigma_guesses = []
        gamma_guesses = []
        amp_guesses = []
        amp_guess_min = []
        amp_guess_max = []
        amp_guess_delta = []
        amp_guesses_range = []

        asym_guesses = []
        asym_guess_delta = []
        sos_guesses = []
        sos_guess_delta = []
        br_guesses = []
        br_guess_delta = []

        #so_guesses = []
        inputPeaks = []
        inputBranching = []
        inputSinglet = []
        inputCK = []
        #inputKE = []
        sigma_guess_delta = []
        gamma_guess_delta = []

        BE_up_lim = []
        BE_low_lim = []
        BE_limit = []
        BE_corr = []
        BE_corr_mult = []

        sigma_up_lim = []
        sigma_low_lim = []
        sigma_limit = []
        sigma_corr = []
        sigma_corr_mult = []

        gamma_up_lim = []
        gamma_low_lim = []
        gamma_limit = []
        gamma_corr = []
        gamma_corr_mult = []

        amp_up_lim = []
        amp_low_lim = []
        amp_limit = []
        amp_corr = []
        amp_corr_mult = []

        asym_up_lim = []
        asym_low_lim = []
        asym_limit = []
        asym_corr = []
        asym_corr_mult = []

        sos_up_lim = []
        sos_low_lim = []
        sos_limit = []
        sos_corr = []
        sos_corr_mult = []

        br_up_lim = []
        br_low_lim = []
        br_limit = []
        br_corr = []
        br_corr_mult = []

       
        #Default values if none are selected in periodic table tab (These have no meaning right now)
        
        global element
        global photoelectronLine
        #Won't throw an error here if added a new peak because the first element is already "N/s"
        try:
            element
        except NameError:
            print("Element Not Selected")
            element =[]
            for i in range(int(self.number_of_peaks.get())):
                element.append('N/s')
        
        try:
            photoelectronLine
        except NameError:
            print("Photoelectron Line Not Selected")
            photoelectronLine = []
            for i in range(int(self.number_of_peaks.get())):
                photoelectronLine.append('N/s')


       
     
        for i in range(int(self.number_of_peaks.get())):
            try:
                elements.append(element[i]) #Added to stop error if user adds a peak after running code 
            except IndexError:
                element.append('N/s')
                elements.append(element[i])
            try:
                photoelectronLines.append(photoelectronLine[i])
            except IndexError:
                photoelectronLine.append('N/s')
                photoelectronLines.append(photoelectronLine[i])
            guesses.append(self.BE_guesses[i].get())
            BE_min = -0.5 #Initially set ranges to +/- 0.5 eV from input value --> This is changed in xps.py
            BE_max = 0.5
            BE_guess_min.append(BE_min)
            BE_guess_max.append(BE_max)
            BE_guess_delta.append(0.01)
            BE_guesses_range.append([BE_min, BE_max, 0.01])
            sigma_guesses.append(self.sigma_guesses[i].get())
            gamma_guesses.append(self.gamma_guesses[i].get())
            amp_guesses.append(self.amp_guesses[i].get())
            amp_min = float(self.amp_guesses[i].get())*-0.10 #AMP range set to +/-10% of input value --> This is larger if baseline is large
            amp_max = float(self.amp_guesses[i].get())*0.10
            amp_guess_min.append(amp_min)
            amp_guess_max.append(amp_max)
            amp_guess_delta.append(0.05)
            amp_guesses_range.append([amp_min, amp_max, 0.05])

            asym_guesses.append(self.asym_guesses[i].get())
            asym_guess_delta.append(0.01)
            sos_guesses.append(self.sos_guesses[i].get())
            sos_guess_delta.append(0.01)
            br_guesses.append(self.br_guesses[i].get())
            br_guess_delta.append(0.01)

            #so_guesses.append(self.so_split[i].get())
            inputPeaks.append(self.peaks[i].get())
            inputBranching.append(self.path_branching[i].get())
            inputSinglet.append(self.peak_singlet[i].get())
            inputCK.append(self.peak_coster_kronig[i].get())
            sigma_guess_delta.append(0.0001)
            gamma_guess_delta.append(0.0001)
          
            BE_up_lim.append(self.BE_up_lim[i].get())
            BE_low_lim.append(self.BE_low_lim[i].get())
            BE_limit.append(self.BE_limit[i].get())
            BE_corr.append(self.BE_corr[i].get())
            BE_corr_mult.append(self.BE_corr_mult[i].get())

            sigma_up_lim.append(self.sigma_up_lim[i].get())
            sigma_low_lim.append(self.sigma_low_lim[i].get())
            sigma_limit.append(self.sigma_limit[i].get())
            sigma_corr.append(self.sigma_corr[i].get())
            sigma_corr_mult.append(self.sigma_corr_mult[i].get())

            gamma_up_lim.append(self.gamma_up_lim[i].get())
            gamma_low_lim.append(self.gamma_low_lim[i].get())
            gamma_limit.append(self.gamma_limit[i].get())
            gamma_corr.append(self.gamma_corr[i].get())
            gamma_corr_mult.append(self.gamma_corr_mult[i].get())

            amp_up_lim.append(self.amp_up_lim[i].get())
            amp_low_lim.append(self.amp_low_lim[i].get())
            amp_limit.append(self.amp_limit[i].get())
            amp_corr.append(self.amp_corr[i].get())
            amp_corr_mult.append(self.amp_corr_mult[i].get())

            asym_up_lim.append(self.asym_up_lim[i].get())
            asym_low_lim.append(self.asym_low_lim[i].get())
            asym_limit.append(self.asym_limit[i].get())
            asym_corr.append(self.asym_corr[i].get())
            asym_corr_mult.append(self.asym_corr_mult[i].get())

            sos_up_lim.append(self.sos_up_lim[i].get())
            sos_low_lim.append(self.sos_low_lim[i].get())
            sos_limit.append(self.sos_limit[i].get())
            sos_corr.append(self.sos_corr[i].get())
            sos_corr_mult.append(self.sos_corr_mult[i].get())

            br_up_lim.append(self.br_up_lim[i].get())
            br_low_lim.append(self.br_low_lim[i].get())
            br_limit.append(self.br_limit[i].get())
            br_corr.append(self.br_corr[i].get())
            br_corr_mult.append(self.br_corr_mult[i].get())

        
        #Expert tab:
        gen_alt = self.gen_late_stage.get()
        alt_lorentz = self.LFWHM_alt
        peak_add_remove = self.data_peak_add

       
           #inputKE.append(self.data_KE[i].get())
        #print("BE UP", BE_up_lim, "BE LOW", BE_low_lim, "BE LIMIT", BE_limit, "BE CORR", BE_corr, )
        #print(amp_guesses_range)



        #\nspinOrbitSplit_range = {spinOrbitSplit_range}
        paths = ("\n\n[Paths]\nnPeaks={nPeaks} \nbackground_type = {bkgn_type} \npeak_type = {peak_type} \nbr_range_min = {br_range_min} \nbr_range_max = {br_range_max} \nbr_range_delta = {br_range_delta} \nbr_guess = {br_guess} \nbr_limited = {br_limited} \nbr_correlated = {br_correlated} \nbr_correlated_mult = {br_correlated_mult} \ngamma_CK_range = {gamma_CK_range} \ngamma_CK_guess = {gamma_CK_guess} \nBE_range_min = {BE_range_min} \nBE_range_max = {BE_range_max} \nBE_range_delta = {BE_range_delta} \nBE = {BE} \nBE_limited = {BE_limited} \nBE_correlated = {BE_correlated} \nBE_correlated_mult = {BE_correlated_mult} \nis_singlet = {is_singlet} \nis_coster_kronig = {is_coster_kronig} \nelement_select = {element_select} \nphotoLine_select = {photoLine_select} \nsos_range_min = {sos_range_min} \nsos_range_max = {sos_range_max} \nsos_range_delta = {sos_range_delta} \nsos_guess = {sos_guess} \nsos_limited = {sos_limited} \nsos_correlated = {sos_correlated} \nsos_correlated_mult = {sos_correlated_mult} \nsigma_range_min = "
                 "{sigma_range_min} \nsigma_range_max = {sigma_range_max} \nsigma_range_delta = {sigma_range_delta} \nsigma_guess = {sigma_guess} \nsigma_limited = {sigma_limited} \nsigma_correlated = {sigma_correlated} \nsigma_correlated_mult = {sigma_correlated_mult} \ngamma_range_min = {gamma_range_min} \ngamma_range_max = {gamma_range_max} \ngamma_range_delta = {gamma_range_delta} \ngamma_guess = {gamma_guess} \ngamma_limited = {gamma_limited} \ngamma_correlated = {gamma_correlated} \ngamma_correlated_mult = {gamma_correlated_mult} \nasym_range_min = {asym_range_min} \nasym_range_max = {asym_range_max} \nasym_range_delta = {asym_range_delta} \nasym_guess = {asym_guess} \nasym_limited = {asym_limited} \nasym_correlated = {asym_correlated} \nasym_correlated_mult = {asym_correlated_mult}  \nasym_CK_range = {asym_CK_range} \nasym_CK_guess = {asym_CK_guess} \nasymmetryDoniach_range = {asymmetryDoniach_range} \namp_range_min = {amp_range_min} \namp_range_max = {amp_range_max} \namp_range_delta = {amp_range_delta} \namp_guess = {amp_guess} \namp_limited = {amp_limited} \namp_correlated = {amp_correlated} \namp_correlated_mult = {amp_correlated_mult} \nk_range = {k_range} \nbackground_range = {background_range} \nslope_range = {slope_range} \ngen_alt_val = {gen_alt_val} \nscale_bool = {scale_bool} \nalt_LFWHM = {alt_LFWHM} \npeak_adding = {peak_adding}"
                 .format(nPeaks=int(self.number_of_peaks.get()),
                         bkgn_type = ",".join(str(i) for i in self.background_types), #Do we need to get this?
                         #bkgn_type = str(self.path_bkgn.get()),
                         peak_type = ",".join(str(i) for i in inputPeaks),
                         #branching_ratio_range=", ".join(str(i) for i in [self.branching_ratio_min, self.branching_ratio_max, self.branching_ratio_delta]),
                         #branching_ratio = ", ".join(str(i) for i in inputBranching),
                         br_range_min =",".join(str(i) for i in br_low_lim),
                         br_range_max =",".join(str(i) for i in br_up_lim),
                         br_range_delta =",".join(str(i) for i in br_guess_delta),
                         br_guess = ", ".join(str(i) for i in br_guesses),
                         br_limited = ", ".join(str(i) for i in br_limit),
                         br_correlated = ", ".join(str(i) for i in br_corr),
                         br_correlated_mult = ", ".join(str(i) for i in br_corr_mult),

                         gamma_CK_range = ", ".join(str(i) for i in [self.gamma_CK_min, self.gamma_CK_max, self.gamma_CK_delta]),
                         gamma_CK_guess = ", ".join(str(i) for i in gamma_guesses), #Making it the same initial inputs as the regular gamma values

                         #BE_range=", ".join(str(i) for i in [self.BE_min, self.BE_max, self.BE_delta]),
                         BE_range_min =",".join(str(i) for i in BE_low_lim),
                         BE_range_max =",".join(str(i) for i in BE_up_lim), #EDIT SIGMA AND GAMMA TO BE LIKE AMP AND BE IN MIN/MAX RANGE
                         BE_range_delta =",".join(str(i) for i in BE_guess_delta),
                         BE = ", ".join(str(i) for i in guesses),
                         BE_limited = ", ".join(str(i) for i in BE_limit),
                         BE_correlated = ", ".join(str(i) for i in BE_corr),
                         BE_correlated_mult = ", ".join(str(i) for i in BE_corr_mult),
                         is_singlet = ", ".join(str(i) for i in inputSinglet),
                         is_coster_kronig = ", ".join(str(i) for i in inputCK),
                         #is_KE = ", ".join(str(i) for i in inputKE),
                         element_select  = ", ".join(str(i) for i in elements),
                         photoLine_select = ", ".join(str(i) for i in photoelectronLines),
                         #spinOrbitSplit_range=", ".join(str(i) for i in [self.spinOrbitSplit_min, self.spinOrbitSplit_max, self.spinOrbitSplit_delta]),
                         #spinOrbitSplit = ", ".join(str(i) for i in so_guesses), #Currently the spin-orbit splitting value is being taken as a constant from the user --> Will have to update to be a part of the GA taking info on each element
                         sos_range_min =",".join(str(i) for i in sos_low_lim),
                         sos_range_max =",".join(str(i) for i in sos_up_lim),
                         sos_range_delta =",".join(str(i) for i in sos_guess_delta),
                         sos_guess = ", ".join(str(i) for i in sos_guesses),
                         sos_limited = ", ".join(str(i) for i in sos_limit),
                         sos_correlated = ", ".join(str(i) for i in sos_corr),
                         sos_correlated_mult = ", ".join(str(i) for i in sos_corr_mult),




                         sigma_range_min =",".join(str(i) for i in sigma_low_lim),
                         sigma_range_max =",".join(str(i) for i in sigma_up_lim),
                         sigma_range_delta =",".join(str(i) for i in sigma_guess_delta),

                         #sigma_range=", ".join(str(i) for i in [self.sigma_min, self.sigma_max, self.sigma_delta]),
                         sigma_guess = ", ".join(str(i) for i in sigma_guesses),
                         sigma_limited = ", ".join(str(i) for i in sigma_limit),
                         sigma_correlated = ", ".join(str(i) for i in sigma_corr),
                         sigma_correlated_mult = ", ".join(str(i) for i in sigma_corr_mult),
                         #fwhm_range=", ".join(str(i) for i in [self.fwhm_min, self.fwhm_max, self.fwhm_delta]),

                         gamma_range_min =",".join(str(i) for i in gamma_low_lim),
                         gamma_range_max =",".join(str(i) for i in gamma_up_lim),
                         gamma_range_delta =",".join(str(i) for i in gamma_guess_delta),
                         

                         gamma_guess = ", ".join(str(i) for i in gamma_guesses),
                         gamma_limited = ", ".join(str(i) for i in gamma_limit),
                         gamma_correlated = ", ".join(str(i) for i in gamma_corr),
                         gamma_correlated_mult = ", ".join(str(i) for i in gamma_corr_mult),

                         amp_range_min =",".join(str(i) for i in amp_low_lim),
                         amp_range_max =",".join(str(i) for i in amp_up_lim),
                         amp_range_delta =",".join(str(i) for i in amp_guess_delta),
                         #amp_range =",".join(str(i) for i in [self.amp_min,self.amp_max,self.amp_delta]),
                         amp_guess = ", ".join(str(i) for i in amp_guesses),
                         amp_limited = ", ".join(str(i) for i in amp_limit),
                         amp_correlated = ", ".join(str(i) for i in amp_corr),
                         amp_correlated_mult = ", ".join(str(i) for i in amp_corr_mult),

                         asym_range_min =",".join(str(i) for i in asym_low_lim),
                         asym_range_max =",".join(str(i) for i in asym_up_lim),
                         asym_range_delta =",".join(str(i) for i in asym_guess_delta),
                         asym_guess = ", ".join(str(i) for i in asym_guesses),
                         asym_limited = ", ".join(str(i) for i in asym_limit),
                         asym_correlated = ", ".join(str(i) for i in asym_corr),
                         asym_correlated_mult = ", ".join(str(i) for i in asym_corr_mult),

                         asym_CK_range = ", ".join(str(i) for i in [self.asym_CK_min, self.asym_CK_max, self.asym_CK_delta]),
                         asym_CK_guess = ", ".join(str(i) for i in asym_guesses), #Making it the same initial inputs as the regular gamma values

                         #asymmetry_range=",".join(str(i) for i in [self.asymmetry_min, self.asymmetry_max, self.asymmetry_delta]),
                         asymmetryDoniach_range=",".join(str(i) for i in [self.asymmetryDoniach_min,self.asymmetryDoniach_max,self.asymmetryDoniach_delta]),
                         k_range = ", ".join(str(i) for i in [0.00,0.1,0.001]), #Does delta have too many decimal points?
                         background_range = ", ".join(str(i) for i in [self.background_min, self.background_max, self.background_delta]),
                         #background_shir_range = ", ".join(str(i) for i in [self.background_shir_min, self.background_shir_max, self.background_shir_delta]),
                         #CTou3_range=", ".join(str(i) for i in [self.C_min, self.C_max, self.C_delta]),
                         #DTou3_range=", ".join(str(i) for i in [self.D_min, self.D_max, self.D_delta]),
                         slope_range = ", ".join(str(i) for i in [self.slope_min, self.slope_max, self.slope_delta]),
                         #exp_amp_range = ", ".join(str(i) for i in [self.exp_amp_min, self.exp_amp_max, self.exp_amp_delta]),
                         #exp_decay_range = ", ".join(str(i) for i in [self.exp_decay_min, self.exp_decay_max, self.exp_decay_delta]),
                         #per_range=", ".join(str(i) for i in [self.percent_min.get(), self.percent_max.get()]),
                         #nu=self.nu.get()
                         gen_alt_val = gen_alt,
                         scale_bool = self.scale_var,
                         alt_LFWHM = alt_lorentz,
                         peak_adding = peak_add_remove
                         ))

        mutations = ("\n\n[Mutations]\nchance_of_mutation = {chance} \noriginal_chance_of_mutation = {original} "
                     "\nmutated_options = {opt}"
                     .format(chance=str(self.chance_of_mutation.get()),
                             original=str(self.original_chance_of_mutation.get()),
                             opt=str(self.mutated_options.get())))
        outputs = ("\n\n[Outputs]\nprint_graph = {graph} \nnum_output_paths = {num} "
                   .format(graph=False, num=False))

        print(str(inputs))
        with open(filename, 'w') as writer:
            writer.write(str(inputs))
            writer.write(str(populations))
            writer.write(str(mutations))
            writer.write(str(paths))
            writer.write(str(outputs))
        return filename

    def loop_gen_ini_same_params(self):
        """
        Will loop through every file in the selected directory and run it with the same parameters
        """
        # file_list = [file.absolute() for file in self.output_folder_path.glob('**/*.ini') if file.is_file()]
        file_list = [filename for filename in self.csv_folder_path.glob('**/*txt') if filename.is_file()]
        # file_list = [filename for filename in self.csv_folder_path.iterdir() if filename.is_file()]
        # stem = self.output_folder_path.stem
        # print(file_list)
        for each in file_list:
            # Change the output file name to match the file name
            name_out = each.stem + '_out.txt'
            # parent = self.output_folder_path.parent
            # print("Output name ", name_out)
            output_path = self.output_folder_path.joinpath(name_out)
            # sets equal so that the name in ini file matches
            self.output_file = output_path
            # Write an ini for for it
            self.csv_generate_from = each
            name = str(each.stem) + '.ini'
            # print("ini name: ", name)
            file_each_path = self.output_folder_path.joinpath(name)
            file_each_path.touch()
            # print('file path: ', str(file_each_path))
            self.write_ini(file_each_path)
        self.multi_known = False

    def generate_directory_popup(self):
        # Popup to ask if the user wants to run all files or just the selected file
        self.multi_known = True
        directory_popup = Toplevel(self.root)
        msg = "Do you want to generate ini files for all files in this directory with the same parameters or just the selected file?"
        entry = ttk.Label(directory_popup, text=msg)
        entry.grid(column=0, row=0, columnspan=2, padx=5, pady=3)
        B1 = ttk.Button(directory_popup, text="All files",
                        command=lambda: [directory_popup.destroy(), self.loop_gen_ini_same_params(),
                                         change_multi_known()])
        B2 = ttk.Button(directory_popup, text="Just this one",
                        command=lambda: [directory_popup.destroy(), self.generate_ini(), change_multi_known()])
        B1.grid(column=0, row=1, padx=5, pady=3, sticky=E)
        B2.grid(column=1, row=1, padx=5, pady=3, sticky=W)
        directory_popup.grid_columnconfigure((0, 1), weight=1)
        directory_popup.grid_rowconfigure((0, 1), weight=1)
        directory_popup.protocol('WM_DELETE_WINDOW', directory_popup.destroy)
        directory_popup.attributes('-topmost', 'true')

        def change_multi_known():
            self.multi_known = False

    def generate_ini(self):
        if self.yes_folder.get() == 1 and self.multi_known is False:  # A folder is selected
            # pop up and ask if generate ini for all files, call appropriate loops
            self.generate_directory_popup()
        else:
            def unique_path():
                counter = 0
                while True:
                    num_name = str(name) + "_" + str(counter) + '.ini'
                    out_name = self.csv_generate_from.stem + "_" + str(counter) + '_out.txt'
                    self.output_file = self.output_folder_path.joinpath(out_name)
                    path = self.output_folder_path.joinpath(num_name)
                    if not path.exists():
                        return path
                    counter += 1

            name = self.csv_generate_from.stem
            file_path = unique_path()
            file_path.touch()
            self.write_ini(file_path)
            return file_path
            # os.chdir(pathlib.Path.cwd().parent)  # change the working directory from gui to nano-indent
            # ini_file = filedialog.asksaveasfilename(initialdir=pathlib.Path.cwd(),
            #                                       title="Choose output ini file",
            #                                      filetypes=[("ini files", "*.ini")])
            # if ini_file is None:
            #   return
            # if isinstance(ini_file, tuple) == False:
            #   if len(ini_file) != 0:
            #      self.write_ini(ini_file)
            #     messagebox.showinfo('', 'Ini file written to {fileloc}'.format(fileloc=ini_file))

            # os.chdir(pathlib.Path.cwd().joinpath('gui'))

    def select_csv_folder(self):
        os.chdir(pathlib.Path.cwd().parent)
        folder_name = pathlib.Path(filedialog.askdirectory(initialdir=pathlib.Path.cwd(), title="Choose a folder"))
        self.csv_folder.set(folder_name)
        self.csv_folder_path = folder_name  # No file has been selected yet - this is the folder
        # This calls a method to create a dropdown menu next to generate ini button of the files in the directory
        self.file_dropdown()
        os.chdir(pathlib.Path.cwd().joinpath('gui'))
        return folder_name

    #def read_input(self, filename):
        # parse with configparser
        # replace the C values in the calibration tab
        #config = configparser.ConfigParser()
        #config.read(filename)
        #self.C0.set(config['Calibrations']['C0'])
        #self.C1.set(config['Calibrations']['C1'])
        #self.C2.set(config['Calibrations']['C2'])
        #self.C3.set(config['Calibrations']['C3'])
        #self.C4.set(config['Calibrations']['C4'])
        #self.C5.set(config['Calibrations']['C5'])
        #self.C6.set(config['Calibrations']['C6'])
        #self.C7.set(config['Calibrations']['C7'])
        #self.C8.set(config['Calibrations']['C8'])
        #self.E_i.set(config['tip_const']['E_i'])
        #self.nu_i.set(config['tip_const']['nu_i'])
        # self.nu.set(config['nu']['nu'])

    #def select_calibration_file(self):
        #os.chdir(pathlib.Path.cwd().parent)  # change the working directory from gui to nano-indent
        #file_name = filedialog.askopenfilename(initialdir=pathlib.Path.cwd(), title="Choose txt/csv", filetypes=(
        #    ("txt files", "*.txt"), ("csv files", "*.csv"), ("all files", "*.*")))
        #self.csv_calibration_file.set(file_name)
        #if file_name:
        #    self.read_input(file_name)
        #os.chdir(pathlib.Path.cwd().joinpath('gui'))
       # return file_name

    def select_output_folder(self):
        """
        Select output folder
        """
        os.chdir(pathlib.Path.cwd().parent)  # change the working directory from gui to nano-indent
        folder_name = pathlib.Path(filedialog.askdirectory(initialdir=pathlib.Path.cwd(), title="Choose a folder"))
        self.output_folder.set(folder_name)
        print("select output folder, folder.get ", self.output_folder.get())
        self.output_folder_path = pathlib.Path(folder_name)
        self.output_file = self.output_folder_path.joinpath('out.txt')
        os.chdir(pathlib.Path.cwd().joinpath('gui'))
        return folder_name

    def loop_direc_same_params(self):
        """
        Will loop through every file in the selected directory and run it with the same parameters
        """
        # file_list = [file.absolute() for file in self.output_folder_path.glob('**/*.ini') if file.is_file()]
        file_list = [filename for filename in self.csv_folder_path.glob('**/*txt') if filename.is_file()]
        # file_list = [filename for filename in self.csv_folder_path.iterdir() if filename.is_file()]
        # stem = self.output_folder_path.stem
        # print(file_list)
        for each in file_list:
            # Change the output file name to match the file name
            name_out = each.stem + '_out.txt'
            # parent = self.output_folder_path.parent
            # print("Output name ", name_out)
            # Create output file
            output_path = self.output_folder_path.joinpath(name_out)
            output_path.touch()
            self.output_file = output_path

            # Write an ini for for it
            self.csv_generate_from = each
            name = str(each.stem) + '.ini'
            print("ini name: ", name)
            file_each_path = self.output_folder_path.joinpath(name)
            file_each_path.touch()
            self.write_ini(file_each_path)
            self.command_list.put(str(file_each_path))
            print("Finished adding to command_list")
            # Run the ini using the run_multi_function
            # self.run_multi_ini()
            # Use labda in button press - next go to run_multi_ini so we can stop

    def loop_direc_diff_params(self):
        # At current this just rins the single file selected
        # Would be cool but a lot of work to make it loop through files and accept new inputs
        # file_list = [filename for filename in self.folder_path.iterdir() if filename.is_file()]
        if not pathlib.Path(self.csv_folder_path.joinpath(self.file_menu.get())).is_file():
            print("No file selected from dropdown menu")
        else:
            name = self.generate_ini()
            self.stop_term()
            command = 'xps_neo -i ' + f'"{name.absolute().as_posix()}"' #Changing to xps_neo
            self.proc = subprocess.Popen("exec " + command, shell=True)

    def directory_popup(self):
        # Popup to ask if the user wants to run all files or just the selected file
        self.multi_known = True
        directory_popup = Toplevel(self.root)
        msg = "Do you want to run all files in this directory with the same parameters or just the selected file?"
        entry = ttk.Label(directory_popup, text=msg)
        entry.grid(column=0, row=0, columnspan=2, padx=5, pady=3)
        B1 = ttk.Button(directory_popup, text="All files",
                        command=lambda: [directory_popup.destroy(), self.loop_direc_same_params(),
                                         self.run_multi_ini()])
        B2 = ttk.Button(directory_popup, text="Just this one",
                        command=lambda: [directory_popup.destroy(), self.loop_direc_diff_params()])
        B1.grid(column=0, row=1, padx=5, pady=3, sticky=E)
        B2.grid(column=1, row=1, padx=5, pady=3, sticky=W)
        directory_popup.grid_columnconfigure((0, 1), weight=1)
        directory_popup.grid_rowconfigure((0, 1), weight=1)
        directory_popup.protocol('WM_DELETE_WINDOW', directory_popup.destroy)
        directory_popup.attributes('-topmost', 'true')

    def run_term(self):
        """
        Runs two separate methods
        if yes folder = 1 means that there is a folder selected
            leads to a popup that allows the user to run all the files with the same parameters
        else a single file is selected and run
        """
        if self.yes_folder.get() == 1:  # A folder is selected
            # pop up and ask if run all files, call appropriate loops
            self.directory_popup()
        else:
            name = self.generate_ini()
            self.stop_term()
            command = 'python ../xps_neo/xps.py -i ' + f'"{name.absolute().as_posix()}"' #changing nano_indent.py to xps.py
            print(command)
            self.proc = subprocess.Popen(''.join(command), shell=True)
            self.proc_list.append(self.proc)

    def build_global(self):
        '''
        Create global tab -  generate ini, run, about, dropdown
        '''

        def about_citation():
            popup = Toplevel()
            popup.wm_title("About: Ver: " + str(self.__version__))
            msg = 'Citation:' \
                  '\nTitle' \
                  '\n Authors' \
                  '\n[Submission], Year'
            cite = ttk.Label(popup, text='Citation:', font='TkTextFont')
            cite.grid(column=0, row=0, sticky=W, padx=self.padx, pady=self.pady)
            citation = scrolledtext.ScrolledText(popup, font="TkTextFont")
            citation.grid(column=0, row=1, padx=self.padx, pady=self.pady)
            with open('media/Citation') as f:
                citation.insert(END, f.read())

            License_Label = ttk.Label(popup, text='License:', font='TkTextFont')
            License_Label.grid(column=0, row=2, sticky=W, padx=self.padx, pady=self.pady)
            license = scrolledtext.ScrolledText(popup)
            license.grid(column=0, row=3, sticky=N + S + W + E, padx=self.padx, pady=self.pady)
            with open('../LICENSE') as f:
                license.insert(END, f.read())
            B1 = ttk.Button(popup, text="Okay", command=popup.destroy)
            B1.grid(column=0, row=4, padx=self.padx, pady=self.pady)

            popup.grid_columnconfigure((1, 3), weight=1)
            popup.grid_rowconfigure((1, 3), weight=1)
            popup.protocol('WM_DELETE_WINDOW', popup.destroy)

        # Column 2 is the dropdown list, is created later as it will only appear if needed
        self.generate_button = ttk.Button(self.root, text="Generate Input", command=self.generate_ini)
        self.generate_button.grid(column=3, row=2, sticky=E, padx=self.padx, pady=self.pady)

        self.run_button = ttk.Button(self.root, text='Run', command=self.run_term)
        self.run_button.grid(column=4, row=2, columnspan=1, sticky=E, padx=self.padx, pady=self.pady)
        self.stop_button = ttk.Button(self.root, text='Stop',
                                      command=lambda: [self.stop_term(), self.run_multi_ini()])
        #self.stop_button = ttk.Button(self.root, text='Stop',command=self.stop_term)
        #How to make XPS Neo actually stop? Is this possible with windows?

        self.stop_button.grid(column=1, row=2, columnspan=1, sticky=W, padx=self.padx, pady=self.pady)

        self.about_button = ttk.Button(self.root, text='About', command=about_citation)
        self.about_button.grid(column=0, row=2, columnspan=1, sticky=W, padx=self.padx, pady=self.pady)

        self.root.grid_columnconfigure(0, weight=1)
        self.root.grid_columnconfigure(1, weight=2)
        self.root.grid_columnconfigure(3, weight=1)

        self.root.grid_rowconfigure(0, weight=1)

        # Create a empty frame
        self.label_frame = LabelFrame(self.root, text="Terminal", padx=5, pady=5)
        self.label_frame.grid(column=0, row=1, columnspan=5, padx=self.padx, pady=self.pady, sticky=E + W + N + S)

        # Create the textbox
        self.label_frame.rowconfigure(0, weight=1)
        self.label_frame.columnconfigure(0, weight=1)
        self.txtbox = scrolledtext.ScrolledText(self.label_frame, width=40, height=5) #Changed height from 10 -Alaina 02/14/2024
        self.txtbox.grid(row=0, column=0, sticky=E + W + N + S)

    def file_dropdown(self):
        p = self.csv_folder_path
        self.filelist = [filename.name for filename in p.glob('**/*txt') if filename.is_file()]
        self.file_menu['values'] = self.filelist
        self.file_menu['state'] = 'readonly'
        # sets the width of the combobox to be the length of the first file in the directory, not perfect but best dynamic solution I could think of
        self.file_menu['width'] = len(self.filelist[0])
        self.file_menu.bind("<<ComboboxSelected>>", self.file_selected)

    def build_inputs_tab(self):
        # Add the tab names
        arr_input = ["Input file", "Input Folder", "Output folder","Lines to skip", "Average Files: ", "Files to Average", "Averaged File", "Load in fit File", "x offset", "y offset"]
        self.description_tabs(arr_input, self.input_tab, row=[3, 4, 5, 6, 8, 9, 10, 13, 14, 15]) #deleted 'calibration file' from row 6

        self.input_tab.grid_columnconfigure(1, weight=1)



        def create_data_obj(x,y,z):
            #separate create func so that it can be called whenever skipLn is updated
            self.expert = self.expert
            self.noise_check = True
            self.range_check = True

            try:
                self.data_obj = data.xps_data(pathlib.Path(self.file_name),int(self.skipLn.get()), float(self.x_offset.get()), float(self.y_offset.get()))
                print("File read successfully")
                #expert_selection = self.checkbutton_expert.get() #How to call button status?
                self.noise_check = self.data_obj.noise_check()
                self.range_check = self.data_obj.range_check()
            except:
                print("Error reading file, try inputting skipped lines")
            
            
            if self.expert == False and self.noise_check == False:
                self.run_button.config(state=DISABLED)
            elif self.expert == False and self.range_check == False:
                self.run_button.config(state=DISABLED)
            else:
                self.run_button.config(state=NORMAL)


            #Put error message for if data is too noisy here

        def select_folder():
            if self.yes_folder.get() == 1:  # When multiple input is checked
                csv_file_button.config(state=DISABLED)
                csv_folder_button.config(state=NORMAL)
                self.csv_folder.set("Please select a folder")
                self.csv_file.set("Folder is selected")
                if self.pertub_check.get() == 1:  # They are also running multiple instances of each file
                    self.checkbutton_whole_folder.config(state='normal')
                # self.file_dropdown()
                # Because no folder is selected yet this errors when put here
            elif self.yes_folder.get() == 0:  # Not Checked
                csv_file_button.config(state=NORMAL)
                csv_folder_button.config(state=DISABLED)
                self.csv_folder.set("File is selected")
                self.csv_file.set("Please select a file")
                if self.pertub_check.get() == 1:  # They are running multiple instances and previously may have selected a folder - disable folder button
                    self.checkbutton_whole_folder.config(state='disabled')
                    self.run_folder.set(
                        False)  # If they previously selected an entire folder need to now only run through the single file

        # functions for input data file
        def select_csv_file():
            os.chdir(pathlib.Path.cwd().parent)  # change the working directory from gui to nano-indent
            file_name = filedialog.askopenfilename(initialdir=pathlib.Path.cwd(), title="Choose txt/csv",
                                                   filetypes=(("txt files", "*.txt"), ("csv files", "*.csv"),
                                                              ("all files", "*.*")))
            if not file_name:
                self.csv_file.set('Please select a file')
            else:
                self.csv_folder.set("File is selected")
                self.csv_file.set(pathlib.Path(file_name))
                self.csv_generate_from = pathlib.Path(file_name)
                self.file_name = file_name
                # create the data objectives

                create_data_obj('','','')

            # disable the dropdown file menu (if user had folder and changed their mind)
            if self.yes_folder.get() == 0:
                self.file_menu.configure(state="disabled")
                self.csv_folder.set("File selected")
            os.chdir(pathlib.Path.cwd().joinpath('gui'))

        def select_csv_folder():
            os.chdir(pathlib.Path.cwd().parent)
            initial_dir = pathlib.Path.cwd()
            folder_name = filedialog.askdirectory(initialdir=pathlib.Path.cwd(), title="Choose a folder")

            if not folder_name:  # They did not select a folder
                self.csv_folder.set("Please choose a folder")
            else:
                folder_path = pathlib.Path(folder_name)
                self.csv_file.set("Folder is selected")
                self.csv_folder.set(folder_path)
                self.csv_folder_path = folder_path  # No file has been selected yet - this is the folder
                # This calls a method to create a dropdown menu next to generate ini button of the files in the directory
                self.file_dropdown()
            os.chdir(pathlib.Path.cwd().joinpath('gui'))

        def get_fit_params():
            #Getting all the values from configuration settings file and updating them in the GUI for running a new fit using the same parameters

            #Can be good for loading in fit files from data published to XPS Oasis. Automate this process? Maybe instead of select fit file, give option to find fit file from online?
            #Want to automate this to also be element/photoelectron line specific. Connect perioidic table to this somehow???


            self.fit_file_selected = True

            


            #fit_file_path = os.path.join(os.getcwd(),args.params)
            config = configparser.ConfigParser()
            print("Reading in fitting parameters from Configuration Setttings file")
            print("Fit file name: ", config.read(self.fit_file.get()))


            #Inputs:
            self.skipLn = config.get('Inputs', 'skipln') #For some reason this still needs to be inputted even if loaded in properly. AttributeError: 'int' object has no attribute 'get_x'
            self.skipLn = IntVar(self.root, self.skipLn)
            self.x_offset = config.get('Inputs', 'x_offset')
            self.x_offset = DoubleVar(self.root, self.x_offset)
            self.y_offset = config.get('Inputs', 'y_offset')
            self.y_offset = DoubleVar(self.root, self.y_offset)



            #Expert tab:
            self.gen_late_stage = config.get('Paths', 'gen_alt_val')
            self.gen_late_stage = IntVar(self.root,self.gen_late_stage)
            
            
            #Reading all the paramters and storing them as StringVar, DoubleVar, etc.
            self.population =   config.get('Populations', 'population')
            self.population = StringVar(self.root, self.population)
            self.num_gen =   config.get('Populations', 'num_gen')
            self.num_gen = StringVar(self.root, self.num_gen)
            self.best_sample =   config.get('Populations', 'best_sample')
            self.best_sample = StringVar(self.root, self.best_sample)
            self.lucky_few =   config.get('Populations', 'lucky_few')
            self.lucky_few = StringVar(self.root, self.lucky_few)
            
            self.chance_of_mutation = config.get('Mutations', 'chance_of_mutation')
            self.chance_of_mutation = StringVar(self.root, self.chance_of_mutation)
            self.original_chance_of_mutation = config.get('Mutations', 'original_chance_of_mutation')
            self.original_chance_of_mutation = StringVar(self.root, self.original_chance_of_mutation)
            self.mutated_options = config.get('Mutations', 'mutated_options')
            self.mutated_options = StringVar(self.root, self.mutated_options)
            
            self.number_of_peaks = config.get('Paths', 'nPeaks')
            self.number_of_peaks = IntVar(self.root,self.number_of_peaks)
            self.background_types = config.get('Paths', 'background_type').split(',')
            peakTypes = config.get('Paths', 'peak_type').split(',')
            #branching_ratio = (config.get('Paths', 'branching_ratio')).split(', ')
            BE_range_min = (config.get('Paths', 'BE_range_min')).split(',')
            BE_range_max = (config.get('Paths', 'BE_range_max')).split(',')
            BE_limited = (config.get('Paths', 'BE_limited')).split(', ')
            BE_correlated = (config.get('Paths', 'BE_correlated')).split(', ')
            BE_correlated_mult = (config.get('Paths', 'BE_correlated_mult')).split(', ') 
            energy = (config.get('Paths', 'BE')).split(', ')
            is_singlet = (config.get('Paths', 'is_singlet')).split(', ')
            is_coster_kronig = (config.get('Paths', 'is_coster_kronig')).split(', ')
            #spinOrbitSplit = (config.get('Paths', 'spinOrbitSplit')).split(', ')

            sigma_range_min = (config.get('Paths', 'sigma_range_min')).split(',')
            sigma_range_max = (config.get('Paths', 'sigma_range_max')).split(',')
            sigma_limited = (config.get('Paths', 'sigma_limited')).split(', ')
            sigma_correlated = (config.get('Paths', 'sigma_correlated')).split(', ')
            sigma_correlated_mult = (config.get('Paths', 'sigma_correlated_mult')).split(', ') 
            sigma_guess = (config.get('Paths', 'sigma_guess')).split(', ')

            gamma_range_min = (config.get('Paths', 'gamma_range_min')).split(',')
            gamma_range_max = (config.get('Paths', 'gamma_range_max')).split(',')
            gamma_limited = (config.get('Paths', 'gamma_limited')).split(', ')
            gamma_correlated = (config.get('Paths', 'gamma_correlated')).split(', ')
            gamma_correlated_mult = (config.get('Paths', 'gamma_correlated_mult')).split(', ') 
            gamma_guess = (config.get('Paths', 'gamma_guess')).split(', ')

            amp_range_min = (config.get('Paths', 'amp_range_min')).split(',')
            amp_range_max = (config.get('Paths', 'amp_range_max')).split(',')
            amp_limited = (config.get('Paths', 'amp_limited')).split(', ')
            amp_correlated = (config.get('Paths', 'amp_correlated')).split(', ')
            amp_correlated_mult = (config.get('Paths', 'amp_correlated_mult')).split(', ') 
            amp_guess = (config.get('Paths', 'amp_guess')).split(', ')

            asym_range_min = (config.get('Paths', 'asym_range_min')).split(',')
            asym_range_max = (config.get('Paths', 'asym_range_max')).split(',')
            asym_limited = (config.get('Paths', 'asym_limited')).split(', ')
            asym_correlated = (config.get('Paths', 'asym_correlated')).split(', ')
            asym_correlated_mult = (config.get('Paths', 'asym_correlated_mult')).split(', ') 
            asym_guess = (config.get('Paths', 'asym_guess')).split(', ')

            sos_range_min = (config.get('Paths', 'sos_range_min')).split(',')
            sos_range_max = (config.get('Paths', 'sos_range_max')).split(',')
            sos_limited = (config.get('Paths', 'sos_limited')).split(', ')
            sos_correlated = (config.get('Paths', 'sos_correlated')).split(', ')
            sos_correlated_mult = (config.get('Paths', 'sos_correlated_mult')).split(', ') 
            sos_guess = (config.get('Paths', 'sos_guess')).split(', ')

            br_range_min = (config.get('Paths', 'br_range_min')).split(',')
            br_range_max = (config.get('Paths', 'br_range_max')).split(',')
            br_limited = (config.get('Paths', 'br_limited')).split(', ')
            br_correlated = (config.get('Paths', 'br_correlated')).split(', ')
            br_correlated_mult = (config.get('Paths', 'br_correlated_mult')).split(', ') 
            br_guess = (config.get('Paths', 'br_guess')).split(', ')
            
            #Updating background checkboxes in Fitting Parameters tab
            for i in self.background_types:
                if i == "SVSC": #changed from Shirley-Sherwood
                    self.shirley_value = IntVar(value=1)
                if i == "Shirley":
                    self.reg_shirley_value = IntVar(value=1)
                if i == "Linear":
                    self.linear_value = IntVar(value=1)
                if i == "Exponential":
                    self.exp_value = IntVar(value=1)
                if i == "2-Param Tougaard":
                    self.toug2_value = IntVar(value=1)
                if i == "3-Param Tougaard":
                    self.toug3_value = IntVar(value=1)
                if i == "Polynomial 1":
                    self.poly1_value = IntVar(value=1)
                if i == "Polynomial 2":
                    self.poly2_value = IntVar(value=1)
                if i == "Polynomial 3":
                    self.poly3_value = IntVar(value=1)
                  
               
            #Converting values into the correct list format for calling in the GA
            for i in range(int(self.number_of_peaks.get())):
                

                #self.peaks[i] = StringVar(self.root,peakTypes[i])
                self.peaks[i].set(peakTypes[i])
                #self.path_branching[i] = DoubleVar(self.root, branching_ratio[i])
                #self.BE_low_lim[i] = DoubleVar(self.root, BE_range_min[i])
                self.BE_low_lim[i].set(BE_range_min[i])
                #self.BE_up_lim[i] = DoubleVar(self.root, BE_range_max[i])
                self.BE_up_lim[i].set(BE_range_max[i])
                #self.BE_limit[i] = BooleanVar(self.root, BE_limited[i])
                self.BE_limit[i].set(BE_limited[i])
                #self.BE_corr[i] = StringVar(self.root, BE_correlated[i])
                self.BE_corr[i].set(BE_correlated[i])
                #self.BE_corr_mult[i] = DoubleVar(self.root, BE_correlated_mult[i])
                self.BE_corr_mult[i].set(BE_correlated_mult[i])
                #self.BE_guesses[i] = DoubleVar(self.root, energy[i])
                self.BE_guesses[i].set(energy[i])
                #self.peak_singlet[i] = BooleanVar(self.root, is_singlet[i])
                self.peak_singlet[i].set(is_singlet[i])
                #self.peak_coster_kronig[i] = BooleanVar(self.root, is_coster_kronig[i])
                self.peak_coster_kronig[i].set(is_coster_kronig[i])
                #self.so_split[i] = DoubleVar(self.root, spinOrbitSplit[i])

                #self.sigma_low_lim[i] = DoubleVar(self.root, sigma_range_min[i])
                self.sigma_low_lim[i].set(sigma_range_min[i])
                #self.sigma_up_lim[i] = DoubleVar(self.root, sigma_range_max[i])
                self.sigma_up_lim[i].set(sigma_range_max[i])
                #self.sigma_limit[i] = BooleanVar(self.root, sigma_limited[i])
                self.sigma_limit[i].set(sigma_limited[i])
                #self.sigma_corr[i] = StringVar(self.root, sigma_correlated[i])
                self.sigma_corr[i].set(sigma_correlated[i])
                #self.sigma_corr_mult[i] = DoubleVar(self.root, sigma_correlated_mult[i])
                self.sigma_corr_mult[i].set(sigma_correlated_mult[i])
                #self.sigma_guesses[i] = DoubleVar(self.root, sigma_guess[i])
                self.sigma_guesses[i].set(sigma_guess[i])

                #self.gamma_low_lim[i] = DoubleVar(self.root, gamma_range_min[i])
                self.gamma_low_lim[i].set(gamma_range_min[i])
                #self.gamma_up_lim[i] = DoubleVar(self.root, gamma_range_max[i])
                self.gamma_up_lim[i].set(gamma_range_max[i])
                #self.gamma_limit[i] = BooleanVar(self.root, gamma_limited[i])
                self.gamma_limit[i].set(gamma_limited[i])
                #self.gamma_corr[i] = StringVar(self.root, gamma_correlated[i])
                self.gamma_corr[i].set(gamma_correlated[i])
                #self.gamma_corr_mult[i] = DoubleVar(self.root, gamma_correlated_mult[i])
                self.gamma_corr_mult[i].set(gamma_correlated_mult[i])
                #self.gamma_guesses[i] = DoubleVar(self.root, gamma_guess[i])
                self.gamma_guesses[i].set(gamma_guess[i])

                #self.amp_low_lim[i] = DoubleVar(self.root, amp_range_min[i])
                self.amp_low_lim[i].set(amp_range_min[i])
                #self.amp_up_lim[i] = DoubleVar(self.root, amp_range_max[i])
                self.amp_up_lim[i].set(amp_range_max[i])
                #self.amp_limit[i] = BooleanVar(self.root, amp_limited[i])
                self.amp_limit[i].set(amp_limited[i])
                #self.amp_corr[i] = StringVar(self.root, amp_correlated[i])
                self.amp_corr[i].set(amp_correlated[i])
                #self.amp_corr_mult[i] = DoubleVar(self.root, amp_correlated_mult[i])
                self.amp_corr_mult[i].set(amp_correlated_mult[i])
                #self.amp_guesses[i] = DoubleVar(self.root, amp_guess[i])
                self.amp_guesses[i].set(amp_guess[i])

                #self.asym_low_lim[i] = DoubleVar(self.root, asym_range_min[i])
                self.asym_low_lim[i].set(asym_range_min[i])
                #self.asym_up_lim[i] = DoubleVar(self.root, asym_range_max[i])
                self.asym_up_lim[i].set(asym_range_max[i])
                #self.asym_limit[i] = BooleanVar(self.root, asym_limited[i])
                self.asym_limit[i].set(asym_limited[i])
                #self.asym_corr[i] = StringVar(self.root, asym_correlated[i])
                self.asym_corr[i].set(asym_correlated[i])
                #self.asym_corr_mult[i] = DoubleVar(self.root, asym_correlated_mult[i])
                self.asym_corr_mult[i].set(asym_correlated_mult[i])
                #self.asym_guesses[i] = DoubleVar(self.root, asym_guess[i])
                self.asym_guesses[i] .set(asym_guess[i])

                #self.sos_low_lim[i] = DoubleVar(self.root, sos_range_min[i])
                self.sos_low_lim[i].set(sos_range_min[i])
                #self.sos_up_lim[i] = DoubleVar(self.root, sos_range_max[i])
                self.sos_up_lim[i].set(sos_range_max[i])
                #self.sos_limit[i] = BooleanVar(self.root, sos_limited[i])
                self.sos_limit[i] = BooleanVar(self.root, sos_limited[i])
                #self.sos_corr[i] = StringVar(self.root, sos_correlated[i])
                self.sos_corr[i].set(sos_correlated[i])
                #self.sos_corr_mult[i] = DoubleVar(self.root, sos_correlated_mult[i])
                self.sos_corr_mult[i].set(sos_correlated_mult[i])
                #self.sos_guesses[i] = DoubleVar(self.root, sos_guess[i])
                self.sos_guesses[i].set(sos_guess[i])

                #self.br_low_lim[i] = DoubleVar(self.root, br_range_min[i])
                self.br_low_lim[i].set(br_range_min[i])
                #self.br_up_lim[i] = DoubleVar(self.root, br_range_max[i])
                self.br_up_lim[i].set(br_range_max[i])
                #self.br_limit[i] = BooleanVar(self.root, br_limited[i])
                self.br_limit[i].set(br_limited[i])
                #self.br_corr[i] = StringVar(self.root, br_correlated[i])
                self.br_corr[i].set(br_correlated[i])
                #self.br_corr_mult[i] = DoubleVar(self.root, br_correlated_mult[i])
                self.br_corr_mult[i].set(br_correlated_mult[i])
                #self.br_guesses[i] = DoubleVar(self.root, br_guess[i])
                self.br_guesses[i].set(br_guess[i])

            #Need to recall these tabs so that the new values are updated in the GUI
            self.build_inputs_tab()
            self.build_population_tab()
            self.build_mutations_tab()
            self.build_fitting_param_tab()
            self.build_param_range_tab()
            self.build_expert_tab()
            create_data_obj('','','')
            

           

            
                
            
           
            


        def select_fit_file():
            
            fit_file_name = filedialog.askopenfilename(initialdir=pathlib.Path.cwd(), title="Choose ini",
                                                   filetypes=(("ini files", "*.ini"),
                                                              ("all files", "*.*")))
            
            if not fit_file_name:
                self.fit_file.set('Please select a file')
            else:
                self.fit_file.set(pathlib.Path(fit_file_name))
                self.fit_generate_from = pathlib.Path(fit_file_name)
                self.fit_file_name = fit_file_name
                get_fit_params() #WHere we get all the values after inmporting the configuration settings file
            
            




        multiple_input_button = ttk.Checkbutton(self.input_tab,
                                                variable=self.yes_folder,
                                                command=select_folder,
                                                offvalue=0, onvalue=1)
        multiple_input_button.grid(column=0, row=2, sticky=E)

        def select_files_to_average():
            os.chdir("..") #change the working directory from gui to XES

            all_files = filedialog.askopenfilenames(initialdir = os.getcwd(), title = "Choose txt/csv", filetypes = (("txt files", "*.txt"),("csv files","*.csv"),("all files","*.*")))

            if not all_files:
                self.all_files.set('Please choose files to average')
                #os.chdir("gui")
            else:
                self.all_files_display.set(all_files)
                self.all_files = all_files
                #os.chdir("gui")
                #print(self.all_files)
            os.chdir(pathlib.Path.cwd().joinpath('gui'))

        def select_file_to_save_averaged_data():
            os.chdir("..")
            averaged_file = filedialog.asksaveasfilename(initialdir = os.getcwd(), title = "Choose txt/csv", filetypes = (("txt files", "*.txt"),("csv files","*.csv"),("all files","*.*")))

            if not averaged_file:
                self.averaged_file.set('Please choose a file to save averaged data')
                # os.chdir("gui")
            else:
                self.averaged_file.set(averaged_file)
                # os.chdir("gui")
            os.chdir(pathlib.Path.cwd().joinpath('gui'))

        def average_selected_data():
            #Finds lines to skip at top of data that are comments
            #lines_to_skip = 2

            # with open(self.all_files[0]) as file:
            #     str=file.readline()
            #     while(not str.__contains__('***')):
            #         lines_to_skip += 1
            #         str=file.readline()
            lines_to_skip = int(self.skipLn.get())

            self.num_files = len(self.all_files)
            self.num_points = len(np.loadtxt(self.all_files[0], skiprows=lines_to_skip, usecols=(0,)))

            self.xPoints_for_data_plot = np.loadtxt(self.all_files[0], skiprows=lines_to_skip, usecols=(0,))
            self.yPoints_for_data_plot = np.loadtxt(self.all_files[0], skiprows=lines_to_skip, usecols=(8,))

            for i in range(self.num_files):
                if i > 0:
                    self.xPoints_for_data_plot += np.array(np.loadtxt(self.all_files[i], skiprows=lines_to_skip, usecols=(0,)))
                    self.yPoints_for_data_plot += np.array(np.loadtxt(self.all_files[i], skiprows=lines_to_skip, usecols=(8,)))

            self.xPoints_avg = np.array(self.xPoints_for_data_plot/self.num_files)
            self.yPoints_avg = np.array(self.yPoints_for_data_plot/self.num_files)

            self.plot_data = np.column_stack((self.xPoints_avg,self.yPoints_avg))
            np.savetxt(self.averaged_file.get(), self.plot_data)

        # Entries:
        # Add the tab entry boxes for inputs
        csv_file_entry = ttk.Entry(self.input_tab, textvariable=self.csv_file, font=self.entryFont)
        csv_file_entry.grid(column=1, row=3, sticky=(W, E))
        csv_folder_entry = ttk.Entry(self.input_tab, textvariable=self.csv_folder, font=self.entryFont)
        csv_folder_entry.grid(column=1, row=4, sticky=(W, E))

        # Add the tab entry boxes for outputs
        output_folder_entry = ttk.Entry(self.input_tab, textvariable=self.output_folder, font=self.entryFont)
        output_folder_entry.grid(column=1, row=5, sticky=(W, E))

        # Add the tab entry boxes for averaging
        entry_files_to_average = ttk.Combobox(self.input_tab, textvariable=self.all_files_display, font=self.entryFont)
        entry_files_to_average.grid(column=1, row=9, sticky=(W,E),padx=self.padx,pady=self.pady)

        entry_averaged_file = ttk.Entry(self.input_tab, textvariable=self.averaged_file, font=self.entryFont)
        entry_averaged_file.grid(column=1, row=10, sticky=(W,E),padx=self.padx,pady=self.pady)

        # Buttons:
        # Adding button to chose file or folder input
        checkbutton_label = ttk.Label(self.input_tab, text="Check to select a folder of input files",
                                      font=self.labelFont)
        checkbutton_label.grid(column=1, row=2, sticky=W)
        
        
        self.KE_check = 2
        def KE_selected():
            global KE_check
            if (self.KE_check % 2) == 0:
                print("Reading data in kinetic energy")
                self.data_KE = True
                
                self.KE_check = 1
            else:
                print("Reading data in binding energy")
                self.data_KE = False

                self.KE_check = 2

        self.checkbutton_KE = ttk.Checkbutton(self.input_tab, text="KE",onvalue= True,offvalue=False, command=KE_selected)
        #self.checkbutton_doublets[i] = ttk.Checkbutton(self.fitting_param_tab, text="Doublet", command=doublet_selected)
        self.checkbutton_KE.grid(column=0, row=2, sticky=W)
        self.checkbutton_KE.state(['!alternate'])


       
        
        self.XES_check = 2
        def XES_selected():
            global XES_check
            if (self.XES_check % 2) == 0:
                print("Reading data in photon energy")
                self.data_XES = True
        
                self.XES_check = 1
            else:
                print("Reading data in binding energy")
                self.data_XES = False

                self.XES_check = 2

        #self.checkbutton_XES = ttk.Checkbutton(self.input_tab, text="XES", onvalue= 0,offvalue=1, command=XES_selected) #Dont need a variable --> Will cause buttons to get mixed up for on/off
        #self.checkbutton_doublets[i] = ttk.Checkbutton(self.fitting_param_tab, text="Doublet", command=doublet_selected)
        #self.checkbutton_XES.grid(column=0, row=12, sticky=W)
        #self.checkbutton_XES.state(['!alternate'])
        
        

        




        # Adding buttons to select each different file/folder
        csv_file_button = ttk.Button(self.input_tab, text="Select File", command=select_csv_file,
                                     style='my.TButton')
        csv_file_button.grid(column=2, row=3, sticky=W)
        # Link this to a variable so when folder is selected file dropdown appears
        csv_folder_button = ttk.Button(self.input_tab, text="Select Folder", command=select_csv_folder,
                                       style='my.TButton')
        csv_folder_button.grid(column=2, row=4, sticky=W)
        csv_folder_button.config(state=DISABLED)  # Unless the multiple file button is checked this will be disabled
        output_folder_button = ttk.Button(self.input_tab, text="Select Folder", command=self.select_output_folder,
                                          style='my.TButton')
        output_folder_button.grid(column=2, row=5, sticky=W)

        # Lines to Skip Entry
        skipLn_entry = ttk.Entry(self.input_tab,textvariable=self.skipLn,font = self.entryFont,width=5)
        skipLn_entry.grid(column=1,row=6,sticky=W)
        self.skipLn.trace_add("write",create_data_obj)

        # X offset Entry
        x_offset_entry = ttk.Entry(self.input_tab,textvariable=self.x_offset,font = self.entryFont,width=5)
        x_offset_entry.grid(column=1,row=14,sticky=W)
        self.x_offset.trace_add("write",create_data_obj)

        # Y offset Entry
        y_offset_entry = ttk.Entry(self.input_tab,textvariable=self.y_offset,font = self.entryFont,width=5)
        y_offset_entry.grid(column=1,row=15,sticky=W)
        self.y_offset.trace_add("write",create_data_obj)

        #Loading in fit file button
        fit_file_button = ttk.Button(self.input_tab, text="Select File", command=select_fit_file,
                                     style='my.TButton')
        fit_file_button.grid(column=2, row=13, sticky=W)

        #Scale Data Button
        self.scale_check = 1
        self.scale = False
        def scale_raw():
           
            if (self.scale_check % 2) == 0: #state of button is off. Used when button is clicked on then off again.
                
                self.scale = False
                self.scale_var = False
                self.scale_check = 1
            else: #state of button is on
                self.scale = True
                self.scale_check = 2
                self.scale_var = True

            print(self.scale)


        #self.checkbutton_scale = ttk.Checkbutton(self.input_tab, text="Scale Data", command=scale_raw)
        #self.checkbutton_scale.grid(column=2, row=6, sticky=W)
        #self.checkbutton_scale.state(['!alternate'])

        #Loading in fit file:
        fit_file_entry = ttk.Entry(self.input_tab, textvariable=self.fit_file, font=self.entryFont)
        fit_file_entry.grid(column=1, row=13, sticky=(W, E))




        # Line between Inputs and Averaging
        separator = ttk.Separator(self.input_tab, orient='horizontal')
        separator.grid(column=0, row=7, columnspan=4, sticky=W + E, padx=self.padx,pady=self.pady)

        # Adding buttons for averaging section
        button_choose_files = ttk.Button(self.input_tab, text="Select Files",command=select_files_to_average, style='my.TButton')
        button_choose_files.grid(column=2, row=9, sticky=W,padx=self.padx,pady=self.pady)

        button_output_file = ttk.Button(self.input_tab,text="Select File",command=select_file_to_save_averaged_data, style='my.TButton')
        button_output_file.grid(column=2, row=10, sticky=W,padx=self.padx,pady=self.pady)

        button_average_data = ttk.Button(self.input_tab,text="Average",command=average_selected_data, style='my.TButton')
        button_average_data.grid(column=0, row=11, sticky=W,padx=self.padx,pady=self.pady)

        #self.file_menu.grid(column=2, row=2, sticky=(W, E))

        #self.file_menu.grid(column=2, row=3, sticky=(W, E))


    def build_population_tab(self):
        """
        Build population tab
        """
        arr_pop = ["Population", "Number of generations", "Best individuals (%)", "Lucky survivor (%)"]
        self.description_tabs(arr_pop, self.population_tab, row=[2, 3, 4, 5])
        population_entry = ttk.Entry(self.population_tab, width=7, textvariable=self.population, font=self.entryFont)
        population_entry.grid(column=2, row=2, sticky=W)
        num_gen_entry = ttk.Entry(self.population_tab, width=7, textvariable=self.num_gen, font=self.entryFont)
        num_gen_entry.grid(column=2, row=3, sticky=W)
        best_sample_entry = ttk.Entry(self.population_tab, width=7, textvariable=self.best_sample, font=self.entryFont)
        best_sample_entry.grid(column=2, row=4, sticky=W)
        lucky_few_entry = ttk.Entry(self.population_tab, width=7, textvariable=self.lucky_few, font=self.entryFont)
        lucky_few_entry.grid(column=2, row=5, sticky=W)

    def build_mutations_tab(self):
        arr_mutations = ["Mutation chance (%)", "Original chance of mutation (%)",
                         "Mutation options"]
        self.description_tabs(arr_mutations, self.mutation_tab, row=[2, 3, 4])
        mut_list = list(range(101))
        chance_of_mutation_entry = ttk.Combobox(self.mutation_tab, width=7, textvariable=self.chance_of_mutation,
                                                values=mut_list,
                                                state="readonly")
        chance_of_mutation_entry.grid(column=4, row=2, sticky=W)
        original_chance_of_mutation_entry = ttk.Combobox(self.mutation_tab, width=7,
                                                         textvariable=self.original_chance_of_mutation,
                                                         values=mut_list, state="readonly")
        original_chance_of_mutation_entry.grid(column=4, row=3, sticky=W)
        mutated_options_drop_list = ttk.Combobox(self.mutation_tab, width=2, textvariable=self.mutated_options,
                                                 values=[0, 1, 2, 3],
                                                 state="readonly")
        mutated_options_drop_list.grid(column=4, row=4, sticky=W)

    
    
    def build_periodicTable_tab(self):
        #Work in progress. Issues with converting from single python file to GUI

        #Wont let me call these in update_element_selection_values so have to give them temp calls here to use
        temp_root = self.root
        BE_guesses = self.BE_guesses
        BE_low_lim = self.BE_low_lim
        gamma_guesses = self.gamma_guesses
        gamma_low_lim = self.gamma_low_lim
        gamma_up_lim = self.gamma_up_lim
        sos_guesses = self.sos_guesses
        br_guesses = self.br_guesses
        singlet = self.peak_singlet
        peakType = self.peaks
        is_coster_kronig = self.peak_coster_kronig
        outer_self = self

        self.photoelectronLineArr = []
        self.elementArr = []
        self.count_table = 0
        for i in range(10):
            self.photoelectronLineArr.append(" ")
            self.elementArr.append(" ")
    
        #build_params = self.build_fitting_param_tab()
        #build_params_range = self.build_param_range_tab()
       

        

      


        #Code taken from https://codereview.stackexchange.com/questions/272438/python-tkinter-periodic-table-of-chemical-elements and adapted to fit XPS data -- Alaina
        #Periodic Table of infomration. When the element is clicked it prompts the user to select the spectral line of interest
        #Used to help narrow down the fitting ranges of parameters in XPS data

        symbols = ['H','He','Li','Be','B','C','N','O','F','Ne',
           'Na','Mg','Al','Si','P','S','Cl','Ar','K', 'Ca',
           'Sc', 'Ti', 'V','Cr', 'Mn', 'Fe', 'Co', 'Ni',
           'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
           'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru',
           'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',
           'I', 'Xe','Cs', 'Ba','La', 'Ce', 'Pr', 'Nd', 'Pm',
           'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm',
           'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir',
           'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
           'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am',
           'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
           'Rf', 'Db', 'Sg', 'Bh','Hs', 'Mt', 'Ds', 'Rg', 'Cn',
           'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']
        keywords =['name','index','element catagory','group','period','block',
                'atomic mass','state of matter','density','electronegativity']

        values = [['Hydrogen',1,'Reactive Nonmetal',1,1,'s',1.01,'gas',0.08,2.2],#H
                 ['Helium',2,'Noble Gas',18,1,'s',4.00,'gas',0.18,'n.A'],#He
                ['Lithium',3,'Alkali Metal',1,2,'s',6.94,'solid',0.53,0.98],#Li
                ['Beryllium',4,'Alkaline Earth Metal',2,2,'s',9.01,'solid',1.84,1.57],#Be
                ['Boron',5,'Metalloid',13,2,'p',10.81,'solid',2.46,2.04],#B
                ['Carbon',6,'Reactive Nonmetal',14,2,'p',12.01,'solid',2.26,2.55],#C
                ['Nitrogen',7,'Reactive Nonmetal',15,2,'p',14.00,'gas',1.17,3.04],#N
                ['Oxygen',8,'Reactive Nonmetal',16,2,'p',15.99,'gas',1.43,3.44],#O
                ['Fluorine',9,'Reactive Nonmetal',17,2,'p',18.99,'gas',1.70,3.98],#F
                ['Neon',10,'Noble Gas',18,2,'p',20.17,'gas',0.90,'n.A'],#Ne

                ['Sodium',11,'Alkali Metal',1,3,'s',22.99,'solid',0.97,0.93],#Na
                ['Magnesium',12,'Alkaline Earth Metal',2,3,'s',24.31,'solid',1.74,1.31],#Mg
                ['Aluminium',13,'Post-transition Metal',13,3,'p',26.98,'solid',2.69,1.61],#Al
                ['Silicon',14,'Metalloid',14,3,'p',28.08,'solid',2.34,1.90],#Si
                ['Phosphorus',15,'Reactive Nonmetal',15,3,'p',30.97,'solid',2.4,2.19],#P
                ['Sulfur',16,'Reactive Nonmetal',16,3,'p',32.06,'solid',2.07,2.58],#S
                ['Chlorine',17,'Reactive Nonmetal',17,3,'p',35.45,'gas',3.22,3.16],#Cl
                ['Argon',18,'Noble Gas',18,3,'p',39.95,'gas',1.78,'n.A'],#Ar
                ['Potassium',19,'Alkali Metal',1,4,'s',39.09,'solid',0.86,0.82],#K
                ['Calicium',20,'Alkaline Earth Metal',2,4,'s',40.08,'solid',1.55,1.00],#Ca

                ['Scandium',21,'Transition Metal',3,4,'d',44.96,'solid',2.99,1.36],#Sc
                ['Titanium',22,'Transition Metal',4,4,'d',47.87,'solid',4.5,1.54],#Ti
                ['Vanadium',23,'Transition Metal',5,4,'d',50.94,'solid',6.11,1.63],#V
                ['Chromium',24,'Transition Metal',6,4,'d',51.99,'solid',7.14,1.66],#Cr
                ['Manganese',25,'Transition Metal',7,4,'d',54.94,'solid',7.43,1.55],#Mn
                ['Iron',26,'Transition Metal',8,4,'d',55.85,'solid',7.87,1.83],#Fe
                ['Cobalt',27,'Transition Metal',9,4,'d',58.93,'solid',8.90,1.88],#Co
                ['Nickel',28,'Transition Metal',10,4,'d',58.69,'solid',8.90,1.91],#Ni

                ['Copper',29,'Transition Metal',11,4,'d',63.54,'solid',8.92,1.90],#Cu
                ['Zinc',30,'Transition Metal',12,4,'d',65.38,'solid',7.14,1.65],#Zn
                ['Gallium',31,'Post-transition Metal',13,4,'p',69.72,'solid',5.90,1.81],#Ga
                ['Germanium',32,'Metalloid',14,4,'p',72.63,'solid',5.32,2.01],#Ge
                ['Arsenic',33,'Metalloid',15,4,'p',74.92,'solid',5.73,2.18],#As
                ['Selenium',34,'Metalloid',16,4,'p',78.97,'solid',4.82,2.55],#Se
                ['Bromine',35,'Reactive Nonmetal',17,4,'p',79.90,'fluid',3.12,2.96],#Br
                ['Krypton',36,'Noble Gas',18,4,'p',83.80,'gas',3.75,3.00],#Kr

                ['Rubidium',37,'Alkali Metal',1,5,'s',85.47,'solid',1.53,0.82],#Rb
                ['Strontium',38,'Alkaline Earth Metal',2,5,'s',87.62,'solid',2.63,0.95],#Sr
                ['Yttrium',39,'Transition Metal',3,5,'d',88.91,'solid',4.47,1.22],#Y
                ['Zirconium',40,'Transition Metal',4,5,'d',91.22,'solid',6.50,1.33],#Zr
                ['Niobium',41,'Transition Metal',5,5,'d',92.90,'solid',8.57,1.6],#Nb
                ['Molybdenum',42,'Transition Metal',6,5,'d',95.95,'solid',10.28,2.16],#Mo
                ['Technetium',43,'Transition Metal',7,5,'d',98.90,'solid',11.5,1.9],#Tc
                ['Ruthenium',44,'Transition Metal',8,5,'d',101.07,'solid',12.37,2.2],#Ru

                ['Rhodium',45,'Transition Metal',9,5,'d',102.90,'solid',12.38,2.28],#Rh
                ['Palladium',46,'Transition Metal',10,5,'d',106.42,'solid',11.99,2.20],#Pd
                ['Silver',47,'Transition Metal',11,5,'d',107.87,'solid',10.49,1.93],#Ag
                ['Cadmium',48,'Transition Metal',12,5,'d',112.41,'solid',8.65,1.69],#Cd
                ['Indium',49,'Post-transition Metal',13,5,'p',114.82,'solid',7.31,1.78],#In
                ['Tin',50,'Post-transition Metal',14,5,'p',118.71,'solid',5.77,1.96],#Sn
                ['Antimony',51,'Metalloid',15,5,'p',121.76,'solid',6.70,2.05],#Sb
                ['Tellurium',52,'Metalloid',16,5,'p',127.60,'solid',6.24,2.10],#Te

                ['Iodine',53,'Reactive Nonmetal',17,5,'p',126.90,'solid',4.94,2.66],#I
                ['Xenon',54,'Noble Gas',18,5,'p',131.29,'gas',5.90,2.6],#Xe
                ['Caesium',55,'Alkali Metal',1,6,'s',132.91,'solid',1.90,0.79],#Cs
                ['Barium',56,'Alkaline Earth Metal',2,6,'s',137.33,'solid',3.62,0.89],#Ba
                ['Lanthanum',57,'Transition Metal',3,6,'d',138.90,'solid',6.17,1.1],#La
                ['Cerium',58,'Lanthanide','La',6,'f',140.12,'solid',6.77,1.12],#Ce
                ['Praseodymium',59,'Lanthanide','La',6,'f',140.91,'solid',6.48,1.13],#Pr
                ['Neodymium',60,'Lanthanide','La',6,'f',144.24,'solid',7.00,1.14],#Nd
                ['Promethium',61,'Lanthanide','La',6,'f',146.91,'solid',7.2,'n.A.'],#Pm

                ['Samarium',62,'Lanthanide','La',6,'f',150.36,'solid',7.54,1.17],#Sm
                ['Europium',63,'Lanthanide','La',6,'f',151.96,'solid',5.25,'n.A'],#Eu
                ['Gadolinium',64,'Lanthanide','La',6,'f',157.25,'solid',7.89,1.20],#Gd
                ['Terbium',65,'Lanthanide','La',6,'f',158.93,'solid',8.25,'n.A'],#Tb
                ['Dysprosium',66,'Lanthanide','La',6,'f',162.50,'solid',8.56,1.22],#Dy
                ['Holmium',67,'Lanthanide','La',6,'f',164.93,'solid',8.78,1.23],#Ho
                ['Erbium',68,'Lanthanide','La',6,'f',167.26,'solid',9.05,1.24],#Er
                ['Thulium',69,'Lanthanide','La',6,'f',168.93,'solid',9.32,1.25],#Tm

                ['Ytterbium',70,'Lanthanide','La',6,'f',173.05,'solid',6.97,'n.A'],#Yb
                ['Lutetium',71,'Lanthanide','La',6,'f',174.97,'solid',9.84,1.27],#Lu
                ['Hafnium',72,'Transition Metal',4,6,'d',178.49,'solid',13.28,1.3],#Hf
                ['Tantalum',73,'Transition Metal',5,6,'d',180.95,'solid',16.65,1.5],#Ta
                ['Tungsten',74,'Transition Metal',6,6,'d',183.84,'solid',19.25,2.36],#W
                ['Rhenium',75,'Transition Metal',7,6,'d',186.21,'solid',21.00,1.9],#Re
                ['Osmium',76,'Transition Metal',8,6,'d',190.23,'solid',22.59,2.2],#Os
                ['Irdium',77,'Transition Metal',9,6,'d',192.22,'solid',22.56,2.2],#Ir

                ['Platinum',78,'Transition Metal',10,6,'d',195.08,'solid',21.45,2.2],#Pt
                ['Gold',79,'Transition Metal',11,6,'d',196.97,'solid',19.32,2.54],#Au
                ['Mercury',80,'Transition Metal',12,6,'d',200.59,'fluid',13.55,2.00],#Hg
                ['Thallium',81,'Post-transition Metal',13,6,'p',204.38,'solid',11.85,1.62],#Tl
                ['Lead',82,'Post-transition Metal',14,6,'p',207.20,'solid',11.34,2.33],#Pb
                ['Bismuth',83,'Post-transition Metal',15,6,'p',208.98,'solid',9.78,2.02],#Bi
                ['Polonium',84,'Post-transition Metal',16,6,'p',209.98,'solid',9.20,2.0],#Po
                ['Astatine',85,'Post-transition Metal',17,6,'p',209.99,'solid','n.A',2.2],#At
                ['Radon',86,'Noble Gas',18,6,'p',222.00,'gas',9.73,'n.A'],#Rn

                ['Francium',87,'Alkali Metal',1,7,'s',223.02,'solid','n.A',0.7],#Fr
                ['Radium',88,'Alkaline Earth Metal',2,7,'s',226.03,'solid',5.5,0.9],#Ra
                ['Actinium',89,'Actinide',3,7,'d',227.03,'solid',10.07,1.1],#Ac
                ['Thorium',90,'Actinide','Ac',7,'f',232.04,'solid',11.72,1.3],#Th
                ['Protactinium',91,'Actinide','Ac',7,'f',231.04,'solid',15.37,1.5],#Pa
                ['Uranium',92,'Actinide','Ac',7,'f',238.03,'solid',19.16,1.38],#U
                ['Neptunium',93,'Actinide','Ac',7,'f',237.05,'solid',20.45,1.36],#Np
                ['Plutonium',94,'Actinide','Ac',7,'f',244.06,'solid',19.82,1.28],#Pu
                ['Americium',95,'Actinide','Ac',7,'f',243.06,'solid',13.67,1.3],#Am

                ['Curium',96,'Actinide','Ac',7,'f',247.07,'solid',13.51,1.3],#Cm
                ['Berkelium',97,'Actinide','Ac',7,'f',247,'solid',14.78,1.3],#Bk
                ['Californium',98,'Actinide','Ac',7,'f',251,'solid',15.1,1.3],#Cf
                ['Einsteinium',99,'Actinide','Ac',7,'f',252,'solid',8.84,'n.A'],#Es
                ['Fermium',100,'Actinide','Ac',7,'f',257.10,'solid','n.A','n.A'],#Fm
                ['Medelevium',101,'Actinide','Ac',7,'f',258,'solid','n.A','n.A'],#Md
                ['Nobelium',102,'Actinide','Ac',7,'f',259,'solid','n.A.','n.A'],#No
                ['Lawrencium',103,'Actinide','Ac',7,'f',266,'solid','n.A','n.A'],#Lr

                ['Rutherfordium',104,'Transition Metal',4,7,'d',261.11,'solid',17.00,'n.A'],#Rf
                ['Dubnium',105,'Transition Metal',5,7,'d',262.11,'n.A','n.A','n.A'],#Db
                ['Seaborgium',106,'Transition Metal',6,7,'d',263.12,'n.A','n.A','n.A'],#Sg
                ['Bohrium',107,'Transition Metal',7,7,'d',262.12,'n.A','n.A','n.A'],#Bh
                ['Hassium',108,'Transition Metal',8,7,'d',265,'n.A','n.A','n.A'],#Hs
                ['Meitnerium',109,'Unknown',9,7,'d',268,'n.A','n.A','n.A'],#Mt
                ['Darmstadtium',110,'Unknown',10,7,'d',281,'n.A','n.A','n.A'],#Ds
                ['Roentgenium',111,'Unknown',11,7,'d',280,'n.A','n.A','n.A'],#Rg
                ['Copernicium',112,'Unknown',12,7,'d',277,'n.A','n.A','n.A'],#Cn

                ['Nihonium',113,'Unknown',13,7,'p',287,'n.A','n.A','n.A'],#Nh
                ['Flerovium',114,'Unknown',14,7,'p',289,'n.A','n.A','n.A'],#Fl
                ['Moscovium',115,'Unknown',15,7,'p',288,'n.A','n.A','n.A'],#Mc
                ['Livermorium',116,'Unknown',16,7,'p',293,'n.A','n.A','n.A'],#Lv
                ['Tennessine',117,'Unknown',17,7,'p',292,'n.A','n.A','n.A'],#Ts
                ['Oganesson',118,'Unknown',18,7,'p',294,'solid',6.6,'n.A']#Og
                ]

        category_colors = {'Alkali Metal' : '#ffabb5',
                            'Alkaline Earth Metal':'#d5b5e6',
                            'Transition Metal':'#91ccff',
                            'Post-transition Metal':'#b6f58c',
                            'Metalloid':'#acc79b',
                            'Reactive Nonmetal':'#f2f18d',
                            'Noble Gas':'#ffc191',
                            'Unknown':'#c8cfca',
                            'Lanthanide':'#a7f3fa',
                            'Actinide':'#a7fade'}
        self.la_offset = -8
        self.ac_offset = -8

        def make_periodicTable(self,symbol,**kwargs):

            
            self.kwargs = kwargs
            self.command= kwargs.pop('command', lambda:print('No command'))
            self.WIDTH,self.HEIGHT,self.BD = 40,40,2
            self.CMP = self.BD*1
            bg = category_colors.get(kwargs.get('element catagory'))


            style = ttk.Style()


            style.configure('table.TFrame', background = bg, foreground = 'black')
            style.configure('table.TLabel', background = bg, foreground = 'black')


            table = ttk.Frame(self.periodicTable_tab, relief = 'raised', style ='table.TFrame' )
            table.configure(width=self.WIDTH,height=self.HEIGHT)
            table.grid_propagate(0)

            self.idx = ttk.Label(table,text=kwargs.get('index'),font=('Arial', 4),style ='table.TLabel') #issue with bg
            #self.u = tk.Label(self,text=kwargs.get('atomic mass'),bg=bg)

            self.name = ttk.Label(table,text=kwargs.get('name'),font=('Arial', 4),style ='table.TLabel')
            symb = ttk.Label(table,text=symbol,font=('bold', 12),style ='table.TLabel')

            #self.e = tk.Label(self,text=kwargs.get('electronegativity'),bg=bg)
            #self.d = tk.Label(self,text=kwargs.get('density'),bg=bg)

            table.grid_columnconfigure(1, weight=2)
            table.grid_rowconfigure(1, weight=2)

            self.idx.grid(row=0,column=0,sticky='w')


            mid_x = self.WIDTH/2-self.name.winfo_reqwidth()/2
            mid_y = self.HEIGHT/2-self.name.winfo_reqheight()/2
            offset= 16
            self.name.place(in_=table,x=mid_x-self.CMP,y=mid_y-self.CMP+offset)

            mid_x = self.WIDTH/2-symb.winfo_reqwidth()/2
            mid_y = self.HEIGHT/2-symb.winfo_reqheight()/2
            symb.place(in_=table,x=mid_x-self.CMP,y=mid_y-self.CMP-offset/2)


            r,c = kwargs.pop('period'),kwargs.pop('group')
            self.offset = 2
            offset = 12
            if c in ('La','Ac'):

                if c == 'La':
                    c =self.la_offset+offset

                    self.la_offset +=1
                    r += self.offset


                if c == 'Ac':
                    c =self.ac_offset+offset
                    self.ac_offset +=1
                    r += offset

            table.grid(row=r,column=c,sticky='nswe')


            
            #Experts tab?
            '''
            self.PT_check = 2
            def PT_selected():
                global PT_check
                if (self.PT_check % 2) == 0:
                    print("Using literature values of selected element")
                    self.data_PT = True
        
                    self.PT_check = 1
                else:
                    print("Using user inputted values")
                    self.data_PT = False

                    self.PT_check = 2

            self.checkbutton_PT = ttk.Checkbutton(self.periodicTable_tab, text="Use Literature Values?", onvalue= 0,offvalue=1, command=PT_selected) #Dont need a variable --> Will cause buttons to get mixed up for on/off
            #self.checkbutton_doublets[i] = ttk.Checkbutton(self.fitting_param_tab, text="Doublet", command=doublet_selected)
            self.checkbutton_PT.grid(column=0, row=35, sticky=W)
            self.checkbutton_PT.state(['!alternate'])
            '''


            def in_active(self):
                #if str(event.type) == 'Enter': self.flag = True
                #if str(event.type) == 'Leave':
                self.flag = False;table.configure(relief='raised')


            def indicate(self): #Want to add in GUI to pop up for spectral line selection

                table.configure(relief='sunken')


            def update_element_selection_values(self):
                #Get values from periodic_table file and update the values in the GUI based on element selection and photoelectron line
                #Updates include: BE value, lorentzian width and width range, if the data should be fit as a doublet and if so the spin-orbit splitting and branching ratio. 
                #If the element selected is a transition metal, peak type defaults to Double Lorentzian, else Voigt selected. Other elements that should use Double Lorentzian?
                #Still under development: Which elements/photoelectron lines should have coster-kronig effects and Tougaard background selection along with other backgrounds

                outer_self.count_table += 1
                global element
                global photoelectronLine
                element = outer_self.elementArr #self.element_select
                photoelectronLine = outer_self.photoelectronLineArr #self.photoLine_select
                self.periodicTable = ElementData(element,photoelectronLine)
                BE_lit, is_singlet, so_split, Br, BE_alt, alt_width, width_range, width, rec_width, default, peakTypes, ck = self.periodicTable.getParams(element,photoelectronLine)

                BE_PT = [0.0]* 10
                width_PT = [0.0]* 10
                width_min = [0.0]* 10
                width_max = [0.0]* 10

              
                for i in range(10):
                
                    if int(round(BE_alt[i])) == 0: 
                        BE_PT[i] = BE_lit[i]
                    else:
                        BE_PT[i] = BE_alt[i]

                    #Default to alt_width, then rec_width, then width
                    if alt_width[i] == 0:
                        width_PT[i] = rec_width[i]
                        if rec_width[i] == 0:
                            width_PT[i] = width[i]
                    else:
                        width_PT[i] = alt_width[i]

                    width_min[i] = -width_range[i]
                    width_max[i] = width_range[i]


                outer_self.background_types = [] #Think it is appending one too many baselines when element is selected. Redefine here to empty list 


                i = 0
                for peak in peakTypes:
                    if peak == "Double Lorentzian":
                        BE_low_lim[i] = (DoubleVar(temp_root, -0.5)) #Allowing peak to go to lower BE because asymmetry makes peak shift to higher BE
                        i +=1
                    

                if outer_self.data_KE == True: #Making sure if data comes in KE that periodic table selection matches KE value
                    for i in range(10):

                        BE_PT[i] = round(1486.6 - BE_PT[i], 2)
         
                
                for i in range(10):
                    
                    #BE_guesses[i] = DoubleVar(temp_root, BE_PT[i]) 
                    BE_guesses[i].set(BE_PT[i])
                    gamma_guesses[i] = DoubleVar(temp_root, width_PT[i]) 
                    gamma_up_lim[i] = DoubleVar(temp_root, width_max[i])
                    gamma_low_lim[i] = DoubleVar(temp_root, width_min[i])
                    br_guesses[i] = DoubleVar(temp_root, Br[i])
                    sos_guesses[i] = DoubleVar(temp_root, so_split[i])
                    singlet[i] = BooleanVar(temp_root, is_singlet[i])
                    peakType[i] = StringVar(temp_root, peakTypes[i])
                    is_coster_kronig[i] = BooleanVar(temp_root, ck[i]) 
             
                
                outer_self.build_fitting_param_tab()
                outer_self.build_param_range_tab()
                
                

            def spectraSelect(self,selection):
                global photoelectronLine
                global element
                element = self.element


                try:
                    outer_self.elementArr[outer_self.count_table] = element
                    outer_self.photoelectronLineArr[outer_self.count_table] = selection
                    print("Peak", outer_self.count_table+1, "assigned:", outer_self.elementArr[outer_self.count_table], outer_self.photoelectronLineArr[outer_self.count_table])
                except:
                    print("Peak maximum reached. Use CLEAR button to reset peak assignment")
                    return
                

                photoelectronLine = selection

                self.element_select = self.element
                self.photoLine_select = selection
                
                update_element_selection_values(self)
            
                '''
                print(self.element, photoelectronLine, "selected")
                if self.element == "H":
                    print("Hydrogen is not detectable using XPS please select a new element")
                if self.element == "He":
                    print("Helium is not detectable in most XPS systems. Make sure you have made the correct selection.")
                '''
                update_peaks_table(0)


            def elementSpectra(self,table):

                top=self.top=Toplevel(table)
                self.element = symb.cget("text")
                photoelectronLine = '1s'
                selectLine = ttk.Label(top, text='Select Photoelectron Line:', font='TkTextFont').pack(side='top',  padx=10,  pady=10)

                if self.element == "H" or self.element == "He" or self.element == "Li" or self.element == "Be" or self.element == "B" or self.element == "C":
                    self.oneS = Button(top,text='1s',command=lambda *args: spectraSelect(self,'1s')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()

                elif self.element == "N" or self.element == "O" or self.element == "F":
                    self.oneS = Button(top,text='1s',command=lambda *args: spectraSelect(self,'1s')).pack(side='left',  padx=10,  pady=10)
                    self.twoS = Button(top,text='2s',command=lambda *args: spectraSelect(self,'2s')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()

                elif self.element == "Ne" or self.element == "Na" or self.element == "Mg":
                    self.oneS = Button(top,text='1s',command=lambda *args: spectraSelect(self,'1s')).pack(side='left',  padx=10,  pady=10)
                    self.twoS = Button(top,text='2s',command=lambda *args: spectraSelect(self,'2s')).pack(side='left',  padx=10,  pady=10)
                    self.twoP = Button(top,text='2p',command=lambda *args: spectraSelect(self,'2p')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()

                elif self.element == "Al" or self.element == "Si":
                    self.twoS = Button(top,text='2s',command=lambda *args: spectraSelect(self,'2s')).pack(side='left',  padx=10,  pady=10)
                    self.twoP = Button(top,text='2p',command=lambda *args: spectraSelect(self,'2p')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()

                elif self.element == "P" or self.element == "S" or self.element == "Ar":
                    self.twoS = Button(top,text='2s',command=lambda *args: spectraSelect(self,'2s')).pack(side='left',  padx=10,  pady=10)
                    self.twoP = Button(top,text='2p',command=lambda *args: spectraSelect(self,'2p')).pack(side='left',  padx=10,  pady=10)
                    self.threeS = Button(top,text='3s',command=lambda *args: spectraSelect(self,'3s')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()

                elif self.element == "Cl" or self.element == "K" or self.element == "Ca" or self.element == "Sc" or self.element == "Ti" or self.element == "V" or self.element == "Cr" or self.element == "Mn" or self.element == "Fe" or self.element == "Co" or self.element == "Ni" or self.element == "Cu":
                    self.twoS = Button(top,text='2s',command=lambda *args: spectraSelect(self,'2s')).pack(side='left',  padx=10,  pady=10)
                    self.twoP = Button(top,text='2p',command=lambda *args: spectraSelect(self,'2p')).pack(side='left',  padx=10,  pady=10)
                    self.threeS = Button(top,text='3s',command=lambda *args: spectraSelect(self,'3s')).pack(side='left',  padx=10,  pady=10)
                    self.threeP = Button(top,text='3p',command=lambda *args: spectraSelect(self,'3p')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()

                elif self.element == "Zn" or self.element == "Ga" or self.element == "Ge" or self.element == "As":
                    self.twoS = Button(top,text='2s',command=lambda *args: spectraSelect(self,'2s')).pack(side='left',  padx=10,  pady=10)
                    self.twoP = Button(top,text='2p',command=lambda *args: spectraSelect(self,'2p')).pack(side='left',  padx=10,  pady=10)
                    self.threeS = Button(top,text='3s',command=lambda *args: spectraSelect(self,'3s')).pack(side='left',  padx=10,  pady=10)
                    self.threeP = Button(top,text='3p',command=lambda *args: spectraSelect(self,'3p')).pack(side='left',  padx=10,  pady=10)
                    self.threeD = Button(top,text='3d',command=lambda *args: spectraSelect(self,'3d')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()

                elif self.element == "Se":
                    self.threeS = Button(top,text='3s',command=lambda *args: spectraSelect(self,'3s')).pack(side='left',  padx=10,  pady=10)
                    self.threeP = Button(top,text='3p',command=lambda *args: spectraSelect(self,'3p')).pack(side='left',  padx=10,  pady=10)
                    self.threeD = Button(top,text='3d',command=lambda *args: spectraSelect(self,'3d')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()
                elif self.element == "Br":
                    self.threeS = Button(top,text='3s',command=lambda *args: spectraSelect(self,'3s')).pack(side='left',  padx=10,  pady=10)
                    self.threeP = Button(top,text='3p',command=lambda *args: spectraSelect(self,'3p')).pack(side='left',  padx=10,  pady=10)
                    self.threeD = Button(top,text='3d',command=lambda *args: spectraSelect(self,'3d')).pack(side='left',  padx=10,  pady=10)
                    self.fourS = Button(top,text='4s',command=lambda *args: spectraSelect(self,'4s')).pack(side='left',  padx=10,  pady=10)
                    self.fourD = Button(top,text='4d',command=lambda *args: spectraSelect(self,'4d')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()

                elif self.element == "Kr" or self.element == "Rb" or self.element == "Sr" or self.element == "Y" or self.element == "Zr" or self.element == "Nb" or self.element == "Mo" or self.element == "Ru" or self.element == "Rh" or self.element == "Pd" or self.element == "Ag":
                    self.threeS = Button(top,text='3s',command=lambda *args: spectraSelect(self,'3s')).pack(side='left',  padx=10,  pady=10)
                    self.threeP = Button(top,text='3p',command=lambda *args: spectraSelect(self,'3p')).pack(side='left',  padx=10,  pady=10)
                    self.threeD = Button(top,text='3d',command=lambda *args: spectraSelect(self,'3d')).pack(side='left',  padx=10,  pady=10)
                    self.fourS = Button(top,text='4s',command=lambda *args: spectraSelect(self,'4s')).pack(side='left',  padx=10,  pady=10)
                    self.fourP = Button(top,text='4p',command=lambda *args: spectraSelect(self,'4p')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()

                elif self.element == "Cd" or self.element == "In" or self.element == "Sb" or self.element == "Te":
                    self.threeS = Button(top,text='3s',command=lambda *args: spectraSelect(self,'3s')).pack(side='left',  padx=10,  pady=10)
                    self.threeP = Button(top,text='3p',command=lambda *args: spectraSelect(self,'3p')).pack(side='left',  padx=10,  pady=10)
                    self.threeD = Button(top,text='3d',command=lambda *args: spectraSelect(self,'3d')).pack(side='left',  padx=10,  pady=10)
                    self.fourS = Button(top,text='4s',command=lambda *args: spectraSelect(self,'4s')).pack(side='left',  padx=10,  pady=10)
                    self.fourP = Button(top,text='4p',command=lambda *args: spectraSelect(self,'4p')).pack(side='left',  padx=10,  pady=10)
                    self.fourD = Button(top,text='4d',command=lambda *args: spectraSelect(self,'4d')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()

                elif self.element == "I" or self.element == "Xe" or self.element == "Cs":
                    self.threeS = Button(top,text='3s',command=lambda *args: spectraSelect(self,'3s')).pack(side='left',  padx=10,  pady=10)
                    self.threeP = Button(top,text='3p',command=lambda *args: spectraSelect(self,'3p')).pack(side='left',  padx=10,  pady=10)
                    self.threeD = Button(top,text='3d',command=lambda *args: spectraSelect(self,'3d')).pack(side='left',  padx=10,  pady=10)
                    self.fourS = Button(top,text='4s',command=lambda *args: spectraSelect(self,'4s')).pack(side='left',  padx=10,  pady=10)
                    self.fourP = Button(top,text='4p',command=lambda *args: spectraSelect(self,'4p')).pack(side='left',  padx=10,  pady=10)
                    self.fourD = Button(top,text='4d',command=lambda *args: spectraSelect(self,'4d')).pack(side='left',  padx=10,  pady=10)
                    self.fiveS = Button(top,text='5s',command=lambda *args: spectraSelect(self,'5s')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()
                elif self.element == "Ba":
                    self.threeS = Button(top,text='3s',command=lambda *args: spectraSelect(self,'3s')).pack(side='left',  padx=10,  pady=10)
                    self.threeP = Button(top,text='3p',command=lambda *args: spectraSelect(self,'3p')).pack(side='left',  padx=10,  pady=10)
                    self.threeD = Button(top,text='3d',command=lambda *args: spectraSelect(self,'3d')).pack(side='left',  padx=10,  pady=10)
                    self.fourS = Button(top,text='4s',command=lambda *args: spectraSelect(self,'4s')).pack(side='left',  padx=10,  pady=10)
                    self.fourP = Button(top,text='4p',command=lambda *args: spectraSelect(self,'4p')).pack(side='left',  padx=10,  pady=10)
                    self.fourD = Button(top,text='4d',command=lambda *args: spectraSelect(self,'4d')).pack(side='left',  padx=10,  pady=10)
                    self.fiveS = Button(top,text='5s',command=lambda *args: spectraSelect(self,'5s')).pack(side='left',  padx=10,  pady=10)
                    self.fiveP = Button(top,text='5p',command=lambda *args: spectraSelect(self,'5p')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()

                elif self.element == "La" or self.element == "Ce" or self.element == "Pr" or self.element == "Nd" or self.element == "Pm":
                    self.threeP = Button(top,text='3p',command=lambda *args: spectraSelect(self,'3p')).pack(side='left',  padx=10,  pady=10)
                    self.threeD = Button(top,text='3d',command=lambda *args: spectraSelect(self,'3d')).pack(side='left',  padx=10,  pady=10)
                    self.fourS = Button(top,text='4s',command=lambda *args: spectraSelect(self,'4s')).pack(side='left',  padx=10,  pady=10)
                    self.fourP = Button(top,text='4p',command=lambda *args: spectraSelect(self,'4p')).pack(side='left',  padx=10,  pady=10)
                    self.fourD = Button(top,text='4d',command=lambda *args: spectraSelect(self,'4d')).pack(side='left',  padx=10,  pady=10)
                    self.fiveS = Button(top,text='5s',command=lambda *args: spectraSelect(self,'5s')).pack(side='left',  padx=10,  pady=10)
                    self.fiveP = Button(top,text='5p',command=lambda *args: spectraSelect(self,'5p')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()

                elif self.element == "Sm" or self.element == "Eu":
                    self.threeD = Button(top,text='3d',command=lambda *args: spectraSelect(self,'3d')).pack(side='left',  padx=10,  pady=10)
                    self.fourS = Button(top,text='4s',command=lambda *args: spectraSelect(self,'4s')).pack(side='left',  padx=10,  pady=10)
                    self.fourP = Button(top,text='4p',command=lambda *args: spectraSelect(self,'4p')).pack(side='left',  padx=10,  pady=10)
                    self.fourD = Button(top,text='4d',command=lambda *args: spectraSelect(self,'4d')).pack(side='left',  padx=10,  pady=10)
                    self.fiveS = Button(top,text='5s',command=lambda *args: spectraSelect(self,'5s')).pack(side='left',  padx=10,  pady=10)
                    self.fiveP = Button(top,text='5p',command=lambda *args: spectraSelect(self,'5p')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()

                elif self.element == "Gd" or self.element == "Tb" or self.element == "Dy" or self.element == "Ho":
                    self.threeD = Button(top,text='3d',command=lambda *args: spectraSelect(self,'3d')).pack(side='left',  padx=10,  pady=10)
                    self.fourS = Button(top,text='4s',command=lambda *args: spectraSelect(self,'4s')).pack(side='left',  padx=10,  pady=10)
                    self.fourP = Button(top,text='4p',command=lambda *args: spectraSelect(self,'4p')).pack(side='left',  padx=10,  pady=10)
                    self.fourD = Button(top,text='4d',command=lambda *args: spectraSelect(self,'4d')).pack(side='left',  padx=10,  pady=10)
                    self.fiveS = Button(top,text='5s',command=lambda *args: spectraSelect(self,'5s')).pack(side='left',  padx=10,  pady=10)
                    self.fiveP = Button(top,text='5p',command=lambda *args: spectraSelect(self,'5p')).pack(side='left',  padx=10,  pady=10)
                    self.fourF = Button(top,text='4f',command=lambda *args: spectraSelect(self,'4f')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()

                elif self.element == "Er" or self.element == "Tm" or self.element == "Yb" or self.element == "Lu" or self.element == "Hf" or self.element == "Ta" or self.element == "W" or self.element == "Os":
                    self.fourS = Button(top,text='4s',command=lambda *args: spectraSelect(self,'4s')).pack(side='left',  padx=10,  pady=10)
                    self.fourP = Button(top,text='4p',command=lambda *args: spectraSelect(self,'4p')).pack(side='left',  padx=10,  pady=10)
                    self.fourD = Button(top,text='4d',command=lambda *args: spectraSelect(self,'4d')).pack(side='left',  padx=10,  pady=10)
                    self.fiveS = Button(top,text='5s',command=lambda *args: spectraSelect(self,'5s')).pack(side='left',  padx=10,  pady=10)
                    self.fiveP = Button(top,text='5p',command=lambda *args: spectraSelect(self,'5p')).pack(side='left',  padx=10,  pady=10)
                    self.fourF = Button(top,text='4f',command=lambda *args: spectraSelect(self,'4f')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()

                elif self.element == "Re" or self.element == "Ir" or self.element == "Pt" or self.element == "Au":
                    self.fourS = Button(top,text='4s',command=lambda *args: spectraSelect(self,'4s')).pack(side='left',  padx=10,  pady=10)
                    self.fourP = Button(top,text='4p',command=lambda *args: spectraSelect(self,'4p')).pack(side='left',  padx=10,  pady=10)
                    self.fourD = Button(top,text='4d',command=lambda *args: spectraSelect(self,'4d')).pack(side='left',  padx=10,  pady=10)
                    self.fiveS = Button(top,text='5s',command=lambda *args: spectraSelect(self,'5s')).pack(side='left',  padx=10,  pady=10)
                    self.fourF = Button(top,text='4f',command=lambda *args: spectraSelect(self,'4f')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()

                elif self.element == "Hg" or self.element == "Tl" or self.element == "Pb" or self.element == "Bi" or self.element == "Po" or self.element == "At" or self.element == "Rn" or self.element == "Fr" or self.element == "Ra" or self.element == "Ac":
                    self.fourS = Button(top,text='4s',command=lambda *args: spectraSelect(self,'4s')).pack(side='left',  padx=10,  pady=10)
                    self.fourP = Button(top,text='4p',command=lambda *args: spectraSelect(self,'4p')).pack(side='left',  padx=10,  pady=10)
                    self.fourD = Button(top,text='4d',command=lambda *args: spectraSelect(self,'4d')).pack(side='left',  padx=10,  pady=10)
                    self.fourF = Button(top,text='4f',command=lambda *args: spectraSelect(self,'4f')).pack(side='left',  padx=10,  pady=10)
                    self.fiveP = Button(top,text='5p',command=lambda *args: spectraSelect(self,'5p')).pack(side='left',  padx=10,  pady=10)
                    self.fiveD = Button(top,text='5d',command=lambda *args: spectraSelect(self,'5d')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()

                elif self.element == "Th" or self.element == "Pa":
                    self.fourS = Button(top,text='4s',command=lambda *args: spectraSelect(self,'4s')).pack(side='left',  padx=10,  pady=10)
                    self.fourP = Button(top,text='4p',command=lambda *args: spectraSelect(self,'4p')).pack(side='left',  padx=10,  pady=10)
                    self.fourD = Button(top,text='4d',command=lambda *args: spectraSelect(self,'4d')).pack(side='left',  padx=10,  pady=10)
                    self.fiveS = Button(top,text='5s',command=lambda *args: spectraSelect(self,'5s')).pack(side='left',  padx=10,  pady=10)
                    self.fourF = Button(top,text='4f',command=lambda *args: spectraSelect(self,'4f')).pack(side='left',  padx=10,  pady=10)
                    self.fiveP = Button(top,text='5p',command=lambda *args: spectraSelect(self,'5p')).pack(side='left',  padx=10,  pady=10)
                    self.fiveD = Button(top,text='5d',command=lambda *args: spectraSelect(self,'5d')).pack(side='left',  padx=10,  pady=10)
                    self.sixS = Button(top,text='6s',command=lambda *args: spectraSelect(self,'6s')).pack(side='left',  padx=10,  pady=10)
                    self.sixP = Button(top,text='6p',command=lambda *args: spectraSelect(self,'6p')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()

                elif self.element == "U":
                    self.fourP = Button(top,text='4p',command=lambda *args: spectraSelect(self,'4p')).pack(side='left',  padx=10,  pady=10)
                    self.fourD = Button(top,text='4d',command=lambda *args: spectraSelect(self,'4d')).pack(side='left',  padx=10,  pady=10)
                    self.fiveS = Button(top,text='5s',command=lambda *args: spectraSelect(self,'5s')).pack(side='left',  padx=10,  pady=10)
                    self.fourF = Button(top,text='4f',command=lambda *args: spectraSelect(self,'4f')).pack(side='left',  padx=10,  pady=10)
                    self.fiveP = Button(top,text='5p',command=lambda *args: spectraSelect(self,'5p')).pack(side='left',  padx=10,  pady=10)
                    self.fiveD = Button(top,text='5d',command=lambda *args: spectraSelect(self,'5d')).pack(side='left',  padx=10,  pady=10)
                    self.sixS = Button(top,text='6s',command=lambda *args: spectraSelect(self,'6s')).pack(side='left',  padx=10,  pady=10)
                    self.sixP = Button(top,text='6p',command=lambda *args: spectraSelect(self,'6p')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()

                elif self.element == "Np":
                    self.fourP = Button(top,text='4p',command=lambda *args: spectraSelect(self,'4p')).pack(side='left',  padx=10,  pady=10)
                    self.fourD = Button(top,text='4d',command=lambda *args: spectraSelect(self,'4d')).pack(side='left',  padx=10,  pady=10)
                    self.fiveS = Button(top,text='5s',command=lambda *args: spectraSelect(self,'5s')).pack(side='left',  padx=10,  pady=10)
                    self.fourF = Button(top,text='4f',command=lambda *args: spectraSelect(self,'4f')).pack(side='left',  padx=10,  pady=10)
                    self.fiveP = Button(top,text='5p',command=lambda *args: spectraSelect(self,'5p')).pack(side='left',  padx=10,  pady=10)
                    self.fiveD = Button(top,text='5d',command=lambda *args: spectraSelect(self,'5d')).pack(side='left',  padx=10,  pady=10)
                    self.sixS = Button(top,text='6s',command=lambda *args: spectraSelect(self,'6s')).pack(side='left',  padx=10,  pady=10)
                    self.sixP = Button(top,text='6p',command=lambda *args: spectraSelect(self,'6p')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()

                elif self.element == "Pu" or self.element == "Am":
                    self.fourP = Button(top,text='4p',command=lambda *args: spectraSelect(self,'4p')).pack(side='left',  padx=10,  pady=10)
                    self.fourD = Button(top,text='4d',command=lambda *args: spectraSelect(self,'4d')).pack(side='left',  padx=10,  pady=10)
                    self.fiveS = Button(top,text='5s',command=lambda *args: spectraSelect(self,'5s')).pack(side='left',  padx=10,  pady=10)
                    self.fourF = Button(top,text='4f',command=lambda *args: spectraSelect(self,'4f')).pack(side='left',  padx=10,  pady=10)
                    self.fiveP = Button(top,text='5p',command=lambda *args: spectraSelect(self,'5p')).pack(side='left',  padx=10,  pady=10)
                    self.fiveD = Button(top,text='5d',command=lambda *args: spectraSelect(self,'5d')).pack(side='left',  padx=10,  pady=10)
                    self.sixS = Button(top,text='6s',command=lambda *args: spectraSelect(self,'6s')).pack(side='left',  padx=10,  pady=10)
                    self.sixP = Button(top,text='6p',command=lambda *args: spectraSelect(self,'6p')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()

                elif self.element == "Cm":
                    self.fourP = Button(top,text='4p',command=lambda *args: spectraSelect(self,'4p')).pack(side='left',  padx=10,  pady=10)
                    self.fiveS = Button(top,text='5s',command=lambda *args: spectraSelect(self,'5s')).pack(side='left',  padx=10,  pady=10)
                    self.fourF = Button(top,text='4f',command=lambda *args: spectraSelect(self,'4f')).pack(side='left',  padx=10,  pady=10)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='left',  padx=5,  pady=5)
                    self.top.wait_window()
                else: #elements Z > 100
                    selectNA = ttk.Label(top, text='No information currently exists in our database for this element', font='TkTextFont').pack(side='top',  padx=5,  pady=5)
                    self.exit = Button(top,text='Exit',command=lambda : top.destroy()).pack(side='top',  padx=5,  pady=5)
                    self.top.wait_window()

                #self.button=Button(top,text='Ok',command=self.cleanup)
                #self.button.pack()
     

            def execute(self):

                #print(symb.cget("text")) #This is the variable we want to pull --> Now we want each each with their respective spectral lines

                elementSpectra(self,table)
                #update_element_selection_values(self)
                #self.b["state"] = "disabled"
                #self.master.wait_window(self.w.top)
                #self.b["state"] = "normal"


                table.configure(relief='raised')


            table.bind('<Enter>', in_active)
            table.bind('<Leave>', in_active)
            table.bind('<ButtonPress-1>', indicate) #Pressing button

            table.bind('<ButtonRelease-1>', execute) #Releasing button
            #print("HELLO")
            [child.bind('<ButtonPress-1>', indicate) for child in table.winfo_children()]
            [child.bind('<ButtonRelease-1>', execute) for child in table.winfo_children()]




        def test():
            print('testing..')



        for idx,symbol in enumerate(symbols):
            kwargs = {}
            for k,v in zip(keywords,values[idx]):
                kwargs.update({k:v})
           
            make_periodicTable(self,symbol,command=test,**kwargs)

        arr_xps_peaks = ["Peak #"]
       
        self.num = int(self.number_of_peaks.get())
        new_num = DoubleVar(self.root, 1)
        element_sel = str(self.element_select[self.num-1].get())

      
        self.description_tabs(arr_xps_peaks, self.periodicTable_tab, row=[30])
        num_peaks_selected = [1,2,3,4,5,6,7,8,9,10]
        arr_xps_table = ["Element", "Assign."]
        self.description_tabs_column(arr_xps_table, self.periodicTable_tab, column=[0, 1], row=31)
     
        def clear_array():
                #Resetting all selections
                global element
                global photoelectronLine
                
                print("Data Cleared")
                outer_self.photoelectronLineArr = [' ']* 10
                outer_self.elementArr = [' ']* 10
    
                element = outer_self.elementArr 
                photoelectronLine = outer_self.photoelectronLineArr 
                
                
                update_peaks_table(0)
                
                self.periodicTable = ElementData(element,photoelectronLine)
                BE_lit, is_singlet, so_split, Br, BE_alt, alt_width, width_range, width, rec_width, default, peakTypes, ck = self.periodicTable.getParams(element,photoelectronLine)

                BE_PT = [0.0]* 10
                width_PT = [0.0]* 10
                width_min = [0.0]* 10
                width_max = [0.0]* 10
                for i in range(10):
                
                    if int(round(BE_alt[i])) == 0: 
                        BE_PT[i] = BE_lit[i]
                    else:
                        BE_PT[i] = BE_alt[i]

                    #Default to alt_width, then rec_width, then width
                    if alt_width[i] == 0:
                        width_PT[i] = rec_width[i]
                        if rec_width[i] == 0:
                            width_PT[i] = width[i]
                    else:
                        width_PT[i] = alt_width[i]

                    width_min[i] = -width_range[i]
                    width_max[i] = width_range[i]
                
                outer_self.background_types = [] #Think it is appending one too many baselines when element is selected. Redefine here to empty list 

                

                for i in range(10):
                    #BE_guesses[i] = DoubleVar(temp_root, BE_PT[i]) 
                    BE_guesses[i].set(BE_PT[i])
                    #gamma_guesses[i] = DoubleVar(temp_root, width_PT[i]) 
                    gamma_guesses[i].set(width_PT[i]) 
                    #gamma_up_lim[i] = DoubleVar(temp_root, width_max[i])
                    gamma_up_lim[i].set(width_max[i])
                    #gamma_low_lim[i] = DoubleVar(temp_root, width_min[i])
                    gamma_low_lim[i].set(width_min[i])
                    #br_guesses[i] = DoubleVar(temp_root, Br[i])
                    br_guesses[i].set(Br[i])
                    #sos_guesses[i] = DoubleVar(temp_root, so_split[i])
                    sos_guesses[i].set(so_split[i])
                    #singlet[i] = BooleanVar(temp_root, is_singlet[i])
                    singlet[i].set(is_singlet[i])
                    #peakType[i] = StringVar(temp_root, peakTypes[i])
                    peakType[i].set(peakTypes[i])
                    #is_coster_kronig[i] = BooleanVar(temp_root, ck[i]) 
                    is_coster_kronig[i].set(ck[i]) 
         
             
                outer_self.build_fitting_param_tab()
                outer_self.build_param_range_tab()
               
                
                outer_self.count_table = 0
               
        def update_peaks_table(args):
            
           
            global element
            global photoelectronLine

            

            self.num = int(new_num.get())
            
            element_sel = str(self.element_select[self.num-1].get())
           

            self.element_entry = ttk.Label(self.periodicTable_tab, text=element[self.num-1], font=self.entryFont)
            self.element_entry.grid(column=0, row=32, sticky=(W, E))
    
            self.line_entry = ttk.Label(self.periodicTable_tab, text=photoelectronLine[self.num-1], font=self.entryFont)
            self.line_entry.grid(column=1, row=32, sticky=(W, E))

        clear_button = ttk.Button(self.periodicTable_tab, text="CLEAR", command=clear_array,style='my.TButton')
        clear_button.grid(column=2, row=30, columnspan=2, sticky=W)
        number_of_peaks_entry_table = ttk.Combobox(self.periodicTable_tab, textvariable=new_num, font=self.entryFont,values= num_peaks_selected, width = 1)
        number_of_peaks_entry_table.grid(column=1, row=30, sticky=(W, E)) #on same row as background checkbox
        number_of_peaks_entry_table.bind('<<ComboboxSelected>>', update_peaks_table)


    def build_fitting_param_tab(self):
        """
        Build fitting parameters tab
        """
       
        arr_xps_peaks = ["No. of Peaks"]
        self.description_tabs(arr_xps_peaks, self.fitting_param_tab, row=[2]) #tabs start at row=1. Row=2 is for background and number_of_peaks


        #Row 3 is for the singlet/doublet checkbuttons

        #Doublet Selected: #Should make this peak dependent in the future

        self.peak_state = 'normal' #initial state for so_split and path_branching




        '''
        checkbutton_doublet = ttk.Checkbutton(self.fitting_param_tab, text="Doublet", command=doublet_selected)
        checkbutton_doublet.grid(column=0, row=3, sticky=W)
        checkbutton_doublet.state(['!alternate']) #initial state of checkbutton is alternate --> need to set it to off otherwise checkbutton is filled with black square
        '''

        #Peak type picker:
        self.peak_types = ['Voigt', 'Gaussian', 'Lorentzian','Double Lorentzian', 'Doniach-Sunjic'] #leave it like this so it can be accessed elsewhere -evan
        path_peakType = ttk.Label(self.fitting_param_tab, text="Fit Type", font=self.labelFont)
        path_peakType.grid_configure(column=2, row=5, sticky=W, padx=self.padx, pady=self.pady)

        BE = ttk.Label(self.fitting_param_tab,text = "Energy",font = self.labelFont)
        BE.grid(column = 3,row = 5,sticky = W, padx= self.padx,pady= self.pady)

        sig = ttk.Label(self.fitting_param_tab,text = "GFWHM",font = self.labelFont)
        sig.grid(column = 4,row = 5,sticky = W, padx= self.padx,pady= self.pady)

        gam = ttk.Label(self.fitting_param_tab,text = "LFWHM",font = self.labelFont)
        gam.grid(column = 5,row = 5,sticky = W, padx= self.padx,pady= self.pady)

        A = ttk.Label(self.fitting_param_tab,text = "Amplitude",font = self.labelFont)
        A.grid(column = 6,row = 5,sticky = W, padx= self.padx,pady= self.pady)


        branching_ratio_button = ttk.Label(self.fitting_param_tab,text = "Area Ratio",font = self.labelFont) #changed from Branching Ratio
        branching_ratio_button.grid(column = 8,row = 5,sticky = W, padx= self.padx,pady= self.pady)

        so_split_button = ttk.Label(self.fitting_param_tab,text = "Spin-Orbit Splitting",font = self.labelFont)
        so_split_button.grid(column = 9,row = 5,sticky = W, padx= self.padx,pady= self.pady)

        self.BE_entries = [0] *10
        self.sigma_range_entries = [0] *10
        self.gamma_range_entries = [0] *10
        self.amp_range_entries = [0] *10
        self.checkbutton_doublets = [0] *10
        self.checkbutton_coster_kronig = [0] *10
        self.peakTypes_entries = [0] *10
        self.branching_entries = [0] *10
        self.so_split_entries = [0] *10 
        self.separators = [0]*10
        self.peak_labels = []

        self.oldNum = 1
        #Creates a row for each peak
      
        def updatePeakSelectionRows(args):
            
            self.num = int(self.number_of_peaks.get())
            peak_labels = []
            rows = []
            i=0

            for row in range(6,6+(2*self.num),2):


                peak_labels.append("Peak " + str(i+1))

                rows.append(row)
                self.val = 1


                self.peakTypes_entries[i] = ttk.Combobox(self.fitting_param_tab, textvariable=self.peaks[i], font=self.entryFont,
                                            values= self.peak_types, width=12)
                self.peakTypes_entries[i].grid(column=2, row=row, sticky=W)

                #Binding Energy:
                self.BE_entries[i] = ttk.Entry(self.fitting_param_tab, textvariable=self.BE_guesses[i], font=self.entryFont, width=4)
                self.BE_entries[i].grid(column=3, row=row, sticky=(W, E))

                self.sigma_range_entries[i] = ttk.Entry(self.fitting_param_tab, textvariable=self.sigma_guesses[i], font=self.entryFont, width=4)
                self.sigma_range_entries[i].grid(column=4, row=row, sticky=(W, E))

                self.gamma_range_entries[i] = ttk.Entry(self.fitting_param_tab, textvariable=self.gamma_guesses[i], font=self.entryFont, width=4)
                self.gamma_range_entries[i].grid(column=5, row=row, sticky=(W, E))

                self.amp_range_entries[i] = ttk.Entry(self.fitting_param_tab, textvariable=self.amp_guesses[i], font=self.entryFont, width=4)
                self.amp_range_entries[i].grid(column=6, row=row, sticky=(W, E))



                self.checkbutton_doublets[i] = ttk.Checkbutton(self.fitting_param_tab, text="Doublet", variable = self.peak_singlet[i],onvalue= False,offvalue=True)
                #self.checkbutton_doublets[i] = ttk.Checkbutton(self.fitting_param_tab, text="Doublet", command=doublet_selected)
                self.checkbutton_doublets[i].grid(column=7, row=row, sticky=W)
                self.checkbutton_doublets[i].state(['!alternate']) #initial state of checkbutton is alternate --> need to set it to off otherwise checkbutton is filled with black square

                self.checkbutton_coster_kronig[i] = ttk.Checkbutton(self.fitting_param_tab, text="Coster-Kronig", variable = self.peak_coster_kronig[i],onvalue= True,offvalue=False)
                self.checkbutton_coster_kronig[i].grid(column=11, row=row, sticky=W)
                self.checkbutton_coster_kronig[i].state(['!alternate'])

                self.branching_entries[i] = ttk.Combobox(self.fitting_param_tab, textvariable=self.br_guesses[i], font=self.entryFont,
                            values=['0.5', '0.666', '0.75'], width=7, state=self.peak_state)
                self.so_split_entries[i] = ttk.Entry(self.fitting_param_tab, textvariable=self.sos_guesses[i], font=self.entryFont, width=4) #added so_split stuff
                self.so_split_entries[i].config(state=self.peak_state)

                #path_branching.config(state='disabled')
                self.branching_entries[i].grid(column=8, row=row, sticky=W)
                self.so_split_entries[i].grid(column=9, row=row, sticky=(W, E))
                self.separators[i] = ttk.Separator(self.fitting_param_tab,orient = 'horizontal')
                self.separators[i].grid(column= 0,row = row+1,columnspan=10,sticky=W+E,padx=self.padx)





               
                i+=1
            #destroy any that were leftover
            for k in np.arange(self.num,self.oldNum):
                #Still having some issues with entry box destruction --> May have something to do with var type of self.num or k range
                # AttributeError: 'int' object has no attribute 'destroy'
                self.peakTypes_entries[k].destroy()
                self.BE_entries[k].destroy()
                self.sigma_range_entries[k].destroy()
                self.gamma_range_entries[k].destroy()
                self.amp_range_entries[k].destroy()
                self.checkbutton_doublets[k].destroy()
                self.checkbutton_coster_kronig[k].destroy()
                self.so_split_entries[k].destroy()
                self.branching_entries[k].destroy()
                self.separators[k].destroy()
                num_diff = self.oldNum - self.num

            self.oldNum = self.num

            for label in self.peak_labels:

                label.destroy()
            #print(peak_labels)

            self.peak_labels = self.description_tabs(peak_labels,self.fitting_param_tab,row = rows)

        #Nunber of Peaks:
        number_of_peaks_options = [1,2,3,4,5,6,7,8,9,10]
        number_of_peaks_entry = ttk.Combobox(self.fitting_param_tab, textvariable=self.number_of_peaks, font=self.entryFont,values= number_of_peaks_options, width = 1)
        number_of_peaks_entry.grid(column=2, row=2, sticky=(W, E)) #on same row as background checkbox
        number_of_peaks_entry.bind('<<ComboboxSelected>>', updatePeakSelectionRows)


        #-------------------------------------Backgrounds-----------------------------------------------------

        #Backgrounds Checkboxes:
        path_bkgn = ttk.Label(self.fitting_param_tab, text="Background Type:", font=self.labelFont)
        path_bkgn.grid_configure(column=3, row=2, sticky=W, padx=self.padx, pady=self.pady)

        global background_types
        #self.background_types = [] #make an array of the background types selected

        #-----Baseline Background------

        self.baseline_selected = 1
        def baseline_bkgn():
            global background_types
            global baseline_selected
            if (self.baseline_selected % 2) == 0: #state of button is off. Used when button is clicked on then off again.
                self.background_types.remove('Baseline')
                self.baseline_selected = 1
            else: #state of button is on
                self.background_types.append('Baseline')
                self.baseline_selected = 2


        defaultOn = IntVar(value=1)
        checkbutton_baseline = ttk.Checkbutton(self.fitting_param_tab,variable = defaultOn, text="Baseline", command=baseline_bkgn)
        checkbutton_baseline.grid(column=5, row=3, sticky=W)
        checkbutton_baseline.state(['!alternate']) #Default is that baseline is always selected

        if self.fit_file_selected == False:  
            baseline_bkgn() #Automatically selecting baseline 
        else:
            for i in self.background_types:
                if i == 'Baseline':
                    self.baseline_selected = 2






     #-----1st Order Polynomial Background------

        self.polynomial1_selected = 1
        def polynomial1_bkgn():
            global background_types
            global polynomial1_selected
            if (self.polynomial1_selected % 2) == 0: #state of button is off. Used when button is clicked on then off again.
                self.background_types.remove('Polynomial 1')
                self.polynomial1_selected = 1
            else: #state of button is on
                self.background_types.append('Polynomial 1')
                self.polynomial1_selected = 2
            
        if self.fit_file_selected == True:  
            for i in self.background_types:
                if i == 'Polynomial 1':
                    self.polynomial1_selected = 2



        checkbutton_poly1 = ttk.Checkbutton(self.fitting_param_tab,variable = self.poly1_value, text="Polynomial 1", command=polynomial1_bkgn)
        checkbutton_poly1.grid(column=4, row=4, sticky=W)
        checkbutton_poly1.state(['!alternate']) #Default is that baseline is always selected



         #-----2nd Order Polynomial Background------

        self.polynomial2_selected = 1
        def polynomial2_bkgn():
            global background_types
            global polynomial2_selected
            if (self.polynomial2_selected % 2) == 0: #state of button is off. Used when button is clicked on then off again.
                self.background_types.remove('Polynomial 2')
                self.polynomial2_selected = 1
            else: #state of button is on
                self.background_types.append('Polynomial 2')
                self.polynomial2_selected = 2

        if self.fit_file_selected == True:  
            for i in self.background_types:
                if i == 'Polynomial 2':
                    self.polynomial2_selected = 2



        checkbutton_poly2 = ttk.Checkbutton(self.fitting_param_tab,variable = self.poly2_value, text="Polynomial 2", command=polynomial2_bkgn)
        checkbutton_poly2.grid(column=5, row=4, sticky=W)
        checkbutton_poly2.state(['!alternate']) #Default is that baseline is always selected




         #-----3nd Order Polynomial Background------

        self.polynomial3_selected = 1
        def polynomial3_bkgn():
            global background_types
            global polynomial3_selected
            if (self.polynomial3_selected % 2) == 0: #state of button is off. Used when button is clicked on then off again.
                self.background_types.remove('Polynomial 3')
                self.polynomial3_selected = 1
            else: #state of button is on
                self.background_types.append('Polynomial 3')
                self.polynomial3_selected = 2

        if self.fit_file_selected == True:  
            for i in self.background_types:
                if i == 'Polynomial 3':
                    self.polynomial3_selected = 2



        checkbutton_poly3 = ttk.Checkbutton(self.fitting_param_tab,variable = self.poly3_value, text="Polynomial 3", command=polynomial3_bkgn)
        checkbutton_poly3.grid(column=6, row=4, sticky=W)
        checkbutton_poly3.state(['!alternate']) #Default is that baseline is always selected


        #-----3 Parameter Tougaard-----

        self.tougaard3_selected = 1
        def tougaard3_bkgn():
            global background_types
            global baseline_selected
            if (self.tougaard3_selected % 2) == 0: #state of button is off. Used when button is clicked on then off again.
                self.background_types.remove('3-Param Tougaard')
                self.tougaard3_selected = 1
            else: #state of button is on
                self.background_types.append('3-Param Tougaard')
                self.tougaard3_selected = 2

        if self.fit_file_selected == True:  
            for i in self.background_types:
                if i == '3-Param Tougaard':
                    self.tougaard3_selected = 2

        checkbutton_tougaard3 = ttk.Checkbutton(self.fitting_param_tab,variable = self.toug3_value,text = "3-Param Tougaard", command = tougaard3_bkgn )
        checkbutton_tougaard3.grid(column=6,row =2,sticky = W)
        checkbutton_tougaard3.state(['!alternate'])





        #-----2 Parameter Tougaard-----

        self.tougaard2_selected = 1
        def tougaard2_bkgn():
            global background_types
            global baseline_selected
            if (self.tougaard2_selected % 2) == 0: #state of button is off. Used when button is clicked on then off again.
                self.background_types.remove('2-Param Tougaard')
                self.tougaard2_selected = 1
            else: #state of button is on
                self.background_types.append('2-Param Tougaard')
                self.tougaard2_selected = 2

        if self.fit_file_selected == True:  
            for i in self.background_types:
                if i == '2-Param Tougaard':
                    self.tougaard2_selected = 2

        checkbutton_tougaard2 = ttk.Checkbutton(self.fitting_param_tab,variable = self.toug2_value,text = "2-Param Tougaard", command = tougaard2_bkgn )
        checkbutton_tougaard2.grid(column=6,row =3,sticky = W)
        checkbutton_tougaard2.state(['!alternate'])





        #-----SVSC Background------ #changed from Shirley-Sherwood

        self.shirley_selected = 1
        def shirley_bkgn(): #This is techincally the SVSC background since it implores different k values for different peaks. Add in regular Shirley?
            global background_types
            global shirley_selected
            if (self.shirley_selected % 2) == 0: #state of button is off. Used when button is clicked on then off again.
                self.background_types.remove('SVSC') #changed from Shirley-Sherwood
                self.shirley_selected = 1
            else: #state of button is on
                self.background_types.append('SVSC') #changed from Shirley-Sherwood
                self.shirley_selected = 2
            return self.shirley_selected
        
        if self.fit_file_selected == True:  
            for i in self.background_types:
                if i == 'SVSC': #changed from Shirley-Sherwood
                    self.shirley_selected = 2
        
        checkbutton_shirley = ttk.Checkbutton(self.fitting_param_tab,variable = self.shirley_value, text="SVSC", command=shirley_bkgn) #changed from Shirley-Sherwood
        checkbutton_shirley.grid(column=4, row=2, sticky=W)
        checkbutton_shirley.state(['!alternate'])


        #Add in regular Shirley 
        self.reg_shirley_selected = 1
        def reg_shirley_bkgn(): #This is techincally the SVSC background since it implores different k values for different peaks. Add in regular Shirley?
            global background_types
            global reg_shirley_selected
            if (self.reg_shirley_selected % 2) == 0: #state of button is off. Used when button is clicked on then off again.
                self.background_types.remove('Shirley') 
                self.reg_shirley_selected = 1
            else: #state of button is on
                self.background_types.append('Shirley') 
                self.reg_shirley_selected = 2
            return self.reg_shirley_selected
        
        if self.fit_file_selected == True:  
            for i in self.background_types:
                if i == 'Shirley':
                    self.reg_shirley_selected = 2
        
        checkbutton_reg_shirley = ttk.Checkbutton(self.fitting_param_tab,variable = self.reg_shirley_value, text="Shirley", command=reg_shirley_bkgn)
        checkbutton_reg_shirley.grid(column=7, row=2, sticky=W)
        checkbutton_reg_shirley.state(['!alternate'])

 

        #-----Integral Slope-----

        self.slope_selected = 1
        def integral_slope_bkgn():
            global background_types
            global slope_selected
            if (self.slope_selected % 2) == 0: #state of button is off. Used when button is clicked on then off again.
                self.background_types.remove('Linear')
                self.slope_selected = 1
            else: #state of button is on
                self.background_types.append('Linear')
                self.slope_selected = 2
            return self.slope_selected
        
        if self.fit_file_selected == True:  
            for i in self.background_types:
                if i == 'Linear':
                    self.slope_selected = 2

        checkbutton_slope = ttk.Checkbutton(self.fitting_param_tab,variable = self.linear_value, text="Linear", command=integral_slope_bkgn)
        checkbutton_slope.grid(column=5, row=2, sticky=W)
        checkbutton_slope.state(['!alternate'])


        #-----Exponential-----
        self.exponential_selected = 1
        def exponential():
            global background_types
            global exponential_selected
            if (self.exponential_selected % 2) == 0: #state of button is off. Used when button is clicked on then off again.
                self.background_types.remove('Exponential')
                self.exponential_selected = 1
            else: #state of button is on
                self.background_types.append('Exponential')
                self.exponential_selected = 2
            return self.exponential_selected
        
        if self.fit_file_selected == True:  
            for i in self.background_types:
                if i == 'Exponential':
                    self.exponential_selected = 2

        checkbutton_exponential = ttk.Checkbutton(self.fitting_param_tab,variable = self.exp_value, text="Exponential", command=exponential)
        checkbutton_exponential.grid(column=4, row=3, sticky=W)
        checkbutton_exponential.state(['!alternate'])





        #Picker for background. No longer using --> changed to checkbuttons to allow for multiple
        #path_bkgns_entry = ttk.Combobox(self.fitting_param_tab, textvariable=self.path_bkgn, font=self.entryFont,
        #                               values=['Shirley-Sherwood', 'Linear']) #need to change to checkButton
        #path_bkgns_entry.grid(column=4, row=2, sticky=W)
        '''
        path_peakType = ttk.Label(self.fitting_param_tab, text="Fit Type", font=self.labelFont)
        path_peakType.grid_configure(column=2, row=5, sticky=W, padx=self.padx, pady=self.pady)
        path_peakTypes_entry = ttk.Combobox(self.fitting_param_tab, textvariable=self.path_peakType, font=self.entryFont,
                                       values= self.peak_types)
        path_peakTypes_entry.grid(column=2, row=5, sticky=W)
        '''



        updatePeakSelectionRows(0)
 













   
    def build_param_range_tab(self):
       
        """
        Build parameter ranges tab
        """
        arr_xps_peaks = ["Peak #"]
       
        self.num = int(self.number_of_peaks.get())
        #self.num = int(self.number_of_peaks.get())
        new_num = DoubleVar(self.root, 1)
        self.description_tabs(arr_xps_peaks, self.param_range_tab, row=[2])
        num_peaks_selected = [1,2,3,4,5,6,7,8,9,10]
        upper_limit_label = ttk.Label(self.param_range_tab,text = "Upper Limit",font = self.labelFont)
        upper_limit_label.grid(column = 3,row = 3,sticky = W, padx= self.padx,pady= self.pady)
        lower_limit_label = ttk.Label(self.param_range_tab,text = "Lower Limit",font = self.labelFont)
        lower_limit_label.grid(column = 2,row = 3,sticky = W, padx= self.padx,pady= self.pady)
        fixed_label = ttk.Label(self.param_range_tab,text = "Fixed",font = self.labelFont)
        fixed_label.grid(column = 4,row = 3,sticky = W, padx= self.padx,pady= self.pady)
        correlated_label = ttk.Label(self.param_range_tab,text = "Correlated Peak #",font = self.labelFont)
        correlated_label.grid(column = 5,row = 3,sticky = W, padx= self.padx,pady= self.pady)
        correlated_mult_label = ttk.Label(self.param_range_tab,text = "Correlated Multiplier",font = self.labelFont)
        correlated_mult_label.grid(column = 6,row = 3,sticky = W, padx= self.padx,pady= self.pady)
        #Maybe change "Correlated" to "Correlated Peak #" and add in "Correlation Multiplier"?

       
        def update_peaks(args):
            #self.num = int(self.number_of_peaks.get()) #This is effecting all the peaks --> Good for assigning correct vars but bad because final value it is left on is considered the number of peaks :(
            self.num = int(new_num.get())
            range_peak = ["Peak " + str(self.num), "Energy", "GFWHM", "LFWHM", "Amplitude", "If Double Lorentzian:", "Asymmetry", "If Doublet:", "Spin-orbit Splitting", "Branching Ratio"]
            
            self.range_peak_labels = self.description_tabs(range_peak,self.param_range_tab,row = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12])

            self.BE_up_lim_entry = ttk.Entry(self.param_range_tab, textvariable=self.BE_up_lim[self.num-1], font=self.entryFont, width=12)
            self.BE_up_lim_entry.grid(column=3, row=4, sticky=(W, E))
            self.BE_low_lim_entry = ttk.Entry(self.param_range_tab, textvariable=self.BE_low_lim[self.num-1], font=self.entryFont, width=12)
            self.BE_low_lim_entry.grid(column=2, row=4, sticky=(W, E))
            self.checkbutton_BE_lim = ttk.Checkbutton(self.param_range_tab, variable = self.BE_limit[self.num-1],onvalue= True,offvalue=False)
            self.checkbutton_BE_lim.grid(column=4, row=4, sticky=(W, E))
            self.checkbutton_BE_lim.state(['!alternate'])
            self.BE_corr_entry = ttk.Combobox(self.param_range_tab, textvariable=self.BE_corr[self.num-1], font=self.entryFont,values= num_peaks_selected, width = 5)
            self.BE_corr_entry.grid(column=5, row=4, sticky=(W, E))
            self.BE_corr_mult_entry = ttk.Entry(self.param_range_tab, textvariable=self.BE_corr_mult[self.num-1], font=self.entryFont, width=12)
            self.BE_corr_mult_entry.grid(column=6, row=4, sticky=(W, E))

            self.sigma_up_lim_entry = ttk.Entry(self.param_range_tab, textvariable=self.sigma_up_lim[self.num-1], font=self.entryFont, width=12)
            self.sigma_up_lim_entry.grid(column=3, row=5, sticky=(W, E))
            self.sigma_low_lim_entry = ttk.Entry(self.param_range_tab, textvariable=self.sigma_low_lim[self.num-1], font=self.entryFont, width=12)
            self.sigma_low_lim_entry.grid(column=2, row=5, sticky=(W, E))
            self.checkbutton_sigma_lim = ttk.Checkbutton(self.param_range_tab, variable = self.sigma_limit[self.num-1],onvalue= True,offvalue=False)
            self.checkbutton_sigma_lim.grid(column=4, row=5, sticky=(W, E))
            self.checkbutton_sigma_lim.state(['!alternate'])
            self.sigma_corr_entry = ttk.Combobox(self.param_range_tab, textvariable=self.sigma_corr[self.num-1], font=self.entryFont,values= num_peaks_selected, width = 5)
            self.sigma_corr_entry.grid(column=5, row=5, sticky=(W, E))
            self.sigma_corr_mult_entry = ttk.Entry(self.param_range_tab, textvariable=self.sigma_corr_mult[self.num-1], font=self.entryFont, width=12)
            self.sigma_corr_mult_entry.grid(column=6, row=5, sticky=(W, E))

            self.gamma_up_lim_entry = ttk.Entry(self.param_range_tab, textvariable=self.gamma_up_lim[self.num-1], font=self.entryFont, width=12)
            self.gamma_up_lim_entry.grid(column=3, row=6, sticky=(W, E))
            self.gamma_low_lim_entry = ttk.Entry(self.param_range_tab, textvariable=self.gamma_low_lim[self.num-1], font=self.entryFont, width=12)
            self.gamma_low_lim_entry.grid(column=2, row=6, sticky=(W, E))
            self.checkbutton_gamma_lim = ttk.Checkbutton(self.param_range_tab, variable = self.gamma_limit[self.num-1],onvalue= True,offvalue=False)
            self.checkbutton_gamma_lim.grid(column=4, row=6, sticky=(W, E))
            self.checkbutton_gamma_lim.state(['!alternate'])
            self.gamma_corr_entry = ttk.Combobox(self.param_range_tab, textvariable=self.gamma_corr[self.num-1], font=self.entryFont,values= num_peaks_selected, width = 5)
            self.gamma_corr_entry.grid(column=5, row=6, sticky=(W, E))
            self.gamma_corr_mult_entry = ttk.Entry(self.param_range_tab, textvariable=self.gamma_corr_mult[self.num-1], font=self.entryFont, width=12)
            self.gamma_corr_mult_entry.grid(column=6, row=6, sticky=(W, E))

            self.amp_up_lim_entry = ttk.Entry(self.param_range_tab, textvariable=self.amp_up_lim[self.num-1], font=self.entryFont, width=12)
            self.amp_up_lim_entry.grid(column=3, row=7, sticky=(W, E))
            self.amp_low_lim_entry = ttk.Entry(self.param_range_tab, textvariable=self.amp_low_lim[self.num-1], font=self.entryFont, width=12)
            self.amp_low_lim_entry.grid(column=2, row=7, sticky=(W, E))
            self.checkbutton_amp_lim = ttk.Checkbutton(self.param_range_tab, variable = self.amp_limit[self.num-1],onvalue= True,offvalue=False)
            self.checkbutton_amp_lim.grid(column=4, row=7, sticky=(W, E))
            self.checkbutton_amp_lim.state(['!alternate'])
            self.amp_corr_entry = ttk.Combobox(self.param_range_tab, textvariable=self.amp_corr[self.num-1], font=self.entryFont,values= num_peaks_selected, width = 5)
            self.amp_corr_entry.grid(column=5, row=7, sticky=(W, E))
            self.amp_corr_mult_entry = ttk.Entry(self.param_range_tab, textvariable=self.amp_corr_mult[self.num-1], font=self.entryFont, width=12)
            self.amp_corr_mult_entry.grid(column=6, row=7, sticky=(W, E))

            self.asym_up_lim_entry = ttk.Entry(self.param_range_tab, textvariable=self.asym_up_lim[self.num-1], font=self.entryFont, width=12)
            self.asym_up_lim_entry.grid(column=3, row=9, sticky=(W, E))
            self.asym_low_lim_entry = ttk.Entry(self.param_range_tab, textvariable=self.asym_low_lim[self.num-1], font=self.entryFont, width=12)
            self.asym_low_lim_entry.grid(column=2, row=9, sticky=(W, E))
            self.checkbutton_asym_lim = ttk.Checkbutton(self.param_range_tab, variable = self.asym_limit[self.num-1],onvalue= True,offvalue=False)
            self.checkbutton_asym_lim.grid(column=4, row=9, sticky=(W, E))
            self.checkbutton_asym_lim.state(['!alternate'])
            self.asym_corr_entry = ttk.Combobox(self.param_range_tab, textvariable=self.asym_corr[self.num-1], font=self.entryFont,values= num_peaks_selected, width = 5)
            self.asym_corr_entry.grid(column=5, row=9, sticky=(W, E))
            self.asym_corr_mult_entry = ttk.Entry(self.param_range_tab, textvariable=self.asym_corr_mult[self.num-1], font=self.entryFont, width=12)
            self.asym_corr_mult_entry.grid(column=6, row=9, sticky=(W, E))

            self.sos_up_lim_entry = ttk.Entry(self.param_range_tab, textvariable=self.sos_up_lim[self.num-1], font=self.entryFont, width=12)
            self.sos_up_lim_entry.grid(column=3, row=11, sticky=(W, E))
            self.sos_low_lim_entry = ttk.Entry(self.param_range_tab, textvariable=self.sos_low_lim[self.num-1], font=self.entryFont, width=12)
            self.sos_low_lim_entry.grid(column=2, row=11, sticky=(W, E))
            self.checkbutton_sos_lim = ttk.Checkbutton(self.param_range_tab, variable = self.sos_limit[self.num-1],onvalue= True,offvalue=False)
            self.checkbutton_sos_lim.grid(column=4, row=11, sticky=(W, E))
            self.checkbutton_sos_lim.state(['!alternate'])
            self.sos_corr_entry = ttk.Combobox(self.param_range_tab, textvariable=self.sos_corr[self.num-1], font=self.entryFont,values= num_peaks_selected, width = 5)
            self.sos_corr_entry.grid(column=5, row=11, sticky=(W, E))
            self.sos_corr_mult_entry = ttk.Entry(self.param_range_tab, textvariable=self.sos_corr_mult[self.num-1], font=self.entryFont, width=12)
            self.sos_corr_mult_entry.grid(column=6, row=11, sticky=(W, E))

            self.br_up_lim_entry = ttk.Entry(self.param_range_tab, textvariable=self.br_up_lim[self.num-1], font=self.entryFont, width=12)
            self.br_up_lim_entry.grid(column=3, row=12, sticky=(W, E))
            self.br_low_lim_entry = ttk.Entry(self.param_range_tab, textvariable=self.br_low_lim[self.num-1], font=self.entryFont, width=12)
            self.br_low_lim_entry.grid(column=2, row=12, sticky=(W, E))
            self.checkbutton_br_lim = ttk.Checkbutton(self.param_range_tab, variable = self.br_limit[self.num-1],onvalue= True,offvalue=False)
            self.checkbutton_br_lim.grid(column=4, row=12, sticky=(W, E))
            self.checkbutton_br_lim.state(['!alternate'])
            self.br_corr_entry = ttk.Combobox(self.param_range_tab, textvariable=self.br_corr[self.num-1], font=self.entryFont,values= num_peaks_selected, width = 5)
            self.br_corr_entry.grid(column=5, row=12, sticky=(W, E))
            self.br_corr_mult_entry = ttk.Entry(self.param_range_tab, textvariable=self.br_corr_mult[self.num-1], font=self.entryFont, width=12)
            self.br_corr_mult_entry.grid(column=6, row=12, sticky=(W, E))
            

          
        #WHY DOES THIS EFFECT THE ARRAY SIZE OF ALL ENTRIES????
        number_of_peaks_entry = ttk.Combobox(self.param_range_tab, textvariable=new_num, font=self.entryFont,values= num_peaks_selected, width = 5)
        number_of_peaks_entry.grid(column=2, row=2, sticky=(W, E)) #on same row as background checkbox
        number_of_peaks_entry.bind('<<ComboboxSelected>>', update_peaks)
        #Need to connect the peak selection Combobox to the BE limit so that it extends to a lower BE if Double Lorentzian
        update_peaks(0)

        #self.range_peak_labels = self.description_tabs(num_peaks_selected,self.param_range_tab,row = [3])
        
        #ADD IN CORRELATED/LIMITED/FIXED RANGES FOR ENERGY, SIGMA, GAMMA, AMPLITUDE














































    def build_plot_tab(self):
        """
        Build plot tab
        """
        self.graph_tab.columnconfigure(0, weight=1)
        self.graph_tab.rowconfigure(1, weight=1)

        def plot_selection():
            #self.data_obj.pre_processing((self.percent_min.get(), self.percent_max.get()))
            data_plot.initial_parameters(self.data_obj, self.data_KE, self.data_XES)
            data_plot.plot_selected()

        def plot_raw():
            global data_KE
            global data_XES
            #self.data_obj.pre_processing((self.percent_min.get(), self.percent_max.get()))
            data_plot.initial_parameters(self.data_obj, self.data_KE, self.data_XES)
            x = []
            y = []
            x.append(self.data_obj.get_x(self.data_KE, self.data_XES))
            y.append(self.data_obj.get_y(self.scale_var))

            data_plot.plot(x,y,'Raw Data', 'Raw Data', self.scale_var,self.data_KE, self.data_XES)
            
        def get_scale():
            return self.scale_var 
           

        def plot_both():
            global data_KE
            global data_XES
            self.data_obj.pre_processing((self.percent_min.get(), self.percent_max.get()))
            data_plot.initial_parameters(self.data_obj, self.data_KE, self.data_XES, title=self.csv_generate_from.stem)
            data_plot.plot_raw_and_selected()

        data_plot = Data_plot(self.graph_tab)
        self.plot_button = ttk.Button(self.graph_tab, text='Plot Data', command=plot_raw)
        self.plot_button.grid(column=0, row=0, columnspan=1, sticky=W, padx=self.padx, pady=self.pady)
        self.plot_selected_button = ttk.Button(self.graph_tab, text='Plot Selected Range', command=plot_selection)
        self.plot_selected_button.grid(column=1, row=0, columnspan=1, sticky=W, padx=self.padx, pady=self.pady)

        self.plot_both_button = ttk.Button(self.graph_tab, text='Plot Raw and Selected', command=plot_both)
        self.plot_both_button.grid(column=2, row=0, columnspan=1, sticky=W, padx=self.padx, pady=self.pady)

        
     
    def generate_randomized_ini(self, multifolder, i):
        pop_range = np.arange(self.pop_min.get(), self.pop_max.get(), 100)
        gen_range = np.arange(self.gen_min.get(), self.gen_max.get(), 5)
        mut_range = np.arange(self.mut_min.get(), self.mut_max.get(), 10)

        self.population.set(np.random.choice(pop_range))
        self.num_gen.set(np.random.choice(gen_range))
        self.chance_of_mutation.set(np.random.choice(mut_range))
        name = self.csv_generate_from.stem

        def unique_path():
            counter = i
            while True:
                num_name = str(name) + "_" + str(counter) + '.ini'
                path = self.output_folder_path.joinpath(num_name)
                if not path.exists():
                    return path
                counter += 1

        file_path = unique_path()
        file_path.touch()
        self.write_ini(file_path)

    def generate_multi_ini(self):
        if self.run_folder.get():  # They want to generate ini for every file in directory
            # generates list of files within the folder
            file_list = [filename for filename in self.csv_folder_path.glob('**/*txt') if filename.is_file()]
            # Loop through each file and generate proper number of ini files
            # stem = self.output_folder_path.stem
            # parent = self.output_folder_path.parent
            for i in range(len(file_list)):
                # set the generate_from file to the current file
                self.csv_generate_from = file_list[i]
                fname = file_list[i].stem
                # create specified number of iterations for this file
                for j in range(self.n_ini.get()):
                    # Gives the output path a unique file name
                    name = fname + '_' + str(j) + '_out' + '.txt'
                    output_path = self.output_folder_path.joinpath(name)
                    # output_path.touch()
                    self.output_file = output_path
                    # print("generate_multi output name: ", self.output_folder_path)
                    self.generate_randomized_ini(self.output_folder_path, j)
        else:
            # stem = self.output_folder_path.stem
            for i in range(self.n_ini.get()):
                # Gives the output path a unique file name
                name = self.csv_generate_from.stem + '_' + str(i) + '_out' + '.txt'
                # parent = self.output_folder_path.parent
                output_path = self.output_folder_path.joinpath(name)
                # output_path.touch()
                self.output_file = output_path
                # print("generate_multi output name: ", self.output_folder_path)
                self.generate_randomized_ini(self.output_folder_path, i)
        return self.output_folder_path

    def stop_all(self):
        self.stop_not_pressed = False
        for i in self.proc_list:
            i.kill()
            i.wait()
        while not self.command_list.empty():
            self.command_list.get()
        print("Stopped xps_neo")

    def run_ini_in_command_list(self, flag):
        global loop
        if self.stop_not_pressed and not self.command_list.empty():
            each = self.command_list.get()
            command = "exec xps_neo -i " + each #changed to xps_neo
            self.proc = subprocess.Popen(command, shell=True)
            self.proc.wait()
            self.proc_list.append(self.proc)
            loop = self.root.after(self.run_ini_in_command_list(True))
        else:
            try:
                self.root.after_cancel(loop)
            except:
                pass

        # while self.stop_not_pressed and len(self.command_list) > 0:

    def set_command_list(self):
        self.output_folder_path = self.generate_multi_ini()
        print("in run multi. flag value", self.stop_not_pressed)
        file_list = [f'"{_file.absolute().as_posix()}"' for _file in self.output_folder_path.glob('**/*.ini') if
                     _file.is_file()]
        print("File list in run multi ", file_list)
        for i in range(len(file_list)):
            print("in for ", i)
            self.command_list.put(file_list[i])

    def run_multi_ini(self):
        # Runs all instances of a single file

        print("before loop")
        global loop
        if self.stop_not_pressed and self.command_list:
            print("Stop not pressed and files exist")
            # print("\n\n\n\n\n")
            each = self.command_list.get()
            command = "exec xps_neo -i " + each #changed from nano_neo to xps_neo
            self.proc = subprocess.Popen(command, shell=True)
            self.proc_list.append(self.proc)
            self.pid_list.append(self.proc.pid)
            # print("Current process ID: ", self.proc.pid)
            loop = self.root.after(0, self.run_multi_ini)
            self.proc.wait()
        else:
            # print("in else")
            try:
                # print("trying to cancel")
                self.root.after_cancel(loop)
                # Empty any unrun commands so they do not run on next iteration
                while self.command_list:
                    self.command_list.get()
            except:
                print("In pass")
                pass
            self.stop_not_pressed = True

        # if self.stop_not_pressed:
        #   print("if self.stop_not_pressed is yes")
        #  self.output_folder_path = self.generate_multi_ini()
        # if self.output_folder.get() =='Please choose a folder to save outputs' or not self.output_folder_path:
        #    print("skipped to not running pls work")
        #   return
        # else:
        #   print("in else of run multi")
        # file_list = [str(filename) for filename in self.output_folder_path.glob('**/*.ini') if filename.is_file()]
        # file_list = [f'"{_file.absolute().as_posix()}"' for _file in self.output_folder_path.glob('**/*.ini') if
        #            _file.is_file()]
        #  print("File list in run multi ", file_list)
        # for i in range(len(file_list)):
        #    print("in for ", i)
        #   self.command_list.append(file_list[i])
        # pls_run(self.stop_not_pressed)
        # else:
        #   print("else of run_multi should cause stop")
        #  pls_run(self.stop_not_pressed) # looks the same but it will send in false (hopefully) & cause after_cancel

        # self.run_ini_in_command_list(True)

    def runningThread(self):
        t1 = Thread(target=self.runningmulti)
        print("In runningThread")
        t1.start()

    def runningmulti(self):
        self.set_command_list()
        self.run_multi_ini()



    def build_expert_tab(self):
        #Tab for changing values in the fitting performance that should only be allowed by an expert
        ready_var = IntVar()


        self.expert_check = 2
        def expert_select():

            if (self.expert_check % 2) == 0:
                print("Expert mode selected")
                self.expert = True
                self.run_button.config(state=NORMAL)
  
                self.expert_check = 1
            else:
                print("Expert mode unselected")
                if self.noise_check == False:
                    self.run_button.config(state=DISABLED)

                elif self.range_check == False:
                    self.run_button.config(state=DISABLED)
                    
                    
                self.expert = False

                self.expert_check = 2
        self.checkbutton_expert = ttk.Checkbutton(self.expert_tab, text="Enable Expert Mode", onvalue= True,offvalue=False, command=expert_select)
        self.checkbutton_expert.grid(column=0, row=0, sticky=(W, E))
        self.checkbutton_expert.state(['!alternate'])

        def select_analysis_folder():
            os.chdir(pathlib.Path.cwd().parent)  # change the working directory from gui to nano-indent

            self.folder_name = filedialog.askdirectory(initialdir=pathlib.Path.cwd(), title="Choose a folder")
            if not self.folder_name:
                self.analysis_dir.set('Please choose a directory')
            else:
                # folder_name = os.path.join(folder_name,'feff')
                self.analysis_dir.set(self.folder_name)

            os.chdir(pathlib.Path.cwd().joinpath('gui'))



        def output_unique():
            params = {
                'base': pathlib.Path.cwd().parent,
                'file': self.csv_generate_from,
                'fileName' : self.csv_file.get(),
                'peaks' : self.peak_types,
                'bkgns' : self.background_types,
                'data obj' : self.data_obj


            }

            self.analysis_obj.unique_parameters(self.analysis_dir,params,self.scale_var,self.data_KE, self.data_XES, peakType='Voigt',title='Fit')


        arr_expert = ["Gen. to alter chi^2", "Alter LFWHM", "Allow Peak Add/Sub"]
        self.description_tabs_column(arr_expert, self.expert_tab, column=[0, 1, 2], row=5)

        #arr_2_expert = ["Output Uniqueness Files"]
        #self.description_tabs_column(arr_2_expert, self.expert_tab, column=[0], row=8)
        
        self.gen_alt_entry = ttk.Entry(self.expert_tab, textvariable=self.gen_late_stage, font=self.entryFont, width=4)
        self.gen_alt_entry.grid(column=0, row=6, sticky=(W, E))

        analysis_folder_button = ttk.Button(self.expert_tab, text="Select Folder for Uniqueness Files",command=select_analysis_folder)  # Add command to export data
        analysis_folder_button.grid(column=0, row=7, sticky=(W, E), padx=self.padx, pady=self.pady, columnspan=2)


        button_plot = ttk.Button(self.expert_tab,text="Output Uniqueness Files",command=output_unique)  # Add command to plot data using postprocessing
        button_plot.grid(column=0, row=8, sticky=(W, E), padx=self.padx, pady=self.pady, columnspan=2)




        self.alt_lorentz_check = 2
        def lorentz_alt_selected():
            global alt_lorentz_check
            if (self.alt_lorentz_check % 2) == 0:
                print("Allowing algorithm to alter LFWHM")
                self.LFWHM_alt = True
        
                self.alt_lorentz_check = 1
            else:
                print("Not allowing algorithm to alter LFWHM")
                self.LFWHM_alt = False

                self.alt_lorentz_check = 2
      
        self.checkbutton_alt_lorentz = ttk.Checkbutton(self.expert_tab, text="Allow alteration", command=lorentz_alt_selected)
        #self.checkbutton_doublets[i] = ttk.Checkbutton(self.fitting_param_tab, text="Doublet", command=doublet_selected)
        self.checkbutton_alt_lorentz.grid(column=1, row=6, sticky=W)
        self.checkbutton_alt_lorentz.state(['!alternate'])





        self.peak_add_check = 2
        def peak_add_selected():
            global peak_add_check
            if (self.peak_add_check % 2) == 0:
                print("Allowing algorithm to add./sub. peaks")
                self.data_peak_add = True
        
                self.peak_add_check = 1
            else:
                print("Not allowing algorithm to add./sub. peaks")
                self.data_peak_add = False

                self.peak_add_check = 2
      
        self.checkbutton_peak_add = ttk.Checkbutton(self.expert_tab, text="Allow peak add./sub.", command=peak_add_selected)
        #self.checkbutton_doublets[i] = ttk.Checkbutton(self.fitting_param_tab, text="Doublet", command=doublet_selected)
        self.checkbutton_peak_add.grid(column=2, row=6, sticky=W)
        self.checkbutton_peak_add.state(['!alternate'])
        
    

    def build_output_tab(self):
        """
        Will allow for multiple iterations over the same data to be performed
        Each time create & save ini, run, save outputs
        """

        # pertub_check = IntVar(self.output_tab, 0)

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
                button_run_nini]
            if self.pertub_check.get() == 0:
                for i in widget_lists:
                    i.config(state='disabled')
                    self.checkbutton_whole_folder.config(state='disabled')
            elif self.pertub_check.get() == 1:
                for i in widget_lists:
                    i.config(state='normal')
                    if self.yes_folder.get() == 1:  # Check to see if folder is selected
                        self.checkbutton_whole_folder.config(state='normal')

        arr_out = ["Print graph", "Steady state exit"]
        self.description_tabs(arr_out, self.output_tab)

        checkbutton_print_graph = ttk.Checkbutton(self.output_tab, var=self.print_graph)
        checkbutton_print_graph.grid(column=1, row=0, sticky=W + E, padx=self.padx)

        checkbutton_steady_state = ttk.Checkbutton(self.output_tab, var=self.steady_state_exit)
        checkbutton_steady_state.grid(column=1, row=1, sticky=W + E, padx=self.padx)

        # Create separators
        separator = ttk.Separator(self.output_tab, orient='horizontal')
        separator.grid(column=0, row=2, columnspan=4, sticky=W + E, padx=self.padx)
        self.output_tab.columnconfigure(3, weight=1)

        arr_out = ["Create Multiple Input Files", "Number of Ini Files", "Pertubations-Population(min,max)",
                   "Pertubations-Generation(min,max)", "Pertubations-Mutation(min,max)"]
        self.description_tabs(arr_out, self.output_tab, row=[3, 5, 6, 7, 8])
        # Create New pertubutuions

        checkbutton_pertub = ttk.Checkbutton(self.output_tab, var=self.pertub_check, command=checkbox_multi)
        checkbutton_pertub.grid(column=1, row=3, sticky=W + E, padx=self.padx)

        pertub_list = list(range(1, 101))

        text = 'Each entry allows user to control perturbation percentage of the desire variables.'
        entry = ttk.Label(self.output_tab, text=text, font=self.labelFont)
        entry.grid_configure(column=0, row=4, columnspan=3, sticky=W + E, padx=self.padx, pady=self.pady)

        entry_n_ini = ttk.Entry(self.output_tab, textvariable=self.n_ini, font=self.entryFont)
        entry_n_ini.grid(column=1, row=5, columnspan=2, sticky=(W, E), padx=self.padx)
        entry_n_ini.config(state='disabled')

        width = 5
        # --------------
        entry_pertub_pop_min = ttk.Entry(self.output_tab, width=width, textvariable=self.pop_min, font=self.entryFont)
        entry_pertub_pop_min.grid(column=1, row=6, sticky=(W, E), padx=self.padx)

        entry_pertub_pop_max = ttk.Entry(self.output_tab, width=width, textvariable=self.pop_max, font=self.entryFont)
        entry_pertub_pop_max.grid(column=2, row=6, sticky=(W, E), padx=self.padx)

        entry_pertub_pop_min.config(state='disabled')
        entry_pertub_pop_max.config(state='disabled')

        # --------------
        entry_pertub_gen_min = ttk.Entry(self.output_tab, width=width, textvariable=self.gen_min, font=self.entryFont)
        entry_pertub_gen_min.grid(column=1, row=7, sticky=(W, E), padx=self.padx)

        entry_pertub_gen_max = ttk.Entry(self.output_tab, width=width, textvariable=self.gen_max, font=self.entryFont)
        entry_pertub_gen_max.grid(column=2, row=7, sticky=(W, E), padx=self.padx)

        entry_pertub_gen_min.config(state='disabled')
        entry_pertub_gen_max.config(state='disabled')

        # --------------
        entry_pertub_mut_min = ttk.Entry(self.output_tab, width=width, textvariable=self.mut_min, font=self.entryFont)
        entry_pertub_mut_min.grid(column=1, row=8, sticky=(W, E), padx=self.padx)

        entry_pertub_mut_max = ttk.Entry(self.output_tab, width=width, textvariable=self.mut_max, font=self.entryFont)
        entry_pertub_mut_max.grid(column=2, row=8, sticky=(W, E), padx=self.padx)

        entry_pertub_mut_min.config(state='disabled')
        entry_pertub_mut_max.config(state='disabled')

        # --------------

        button_gen_nini = ttk.Button(self.output_tab, text="Generate Input Files", command=self.generate_multi_ini)
        button_gen_nini.grid(column=0, row=9, columnspan=3, sticky=W + E, padx=self.padx, pady=self.pady)
        button_gen_nini.config(state='disabled')

        button_run_nini = ttk.Button(self.output_tab, text="Run All Instances",
                                     command=lambda: [self.set_command_list(), self.run_multi_ini()])
        button_run_nini.grid(column=0, row=10, columnspan=3, sticky=W + E, padx=self.padx, pady=self.pady)
        button_run_nini.config(state='disabled')

        # Adding button to chose if run all files in the folder
        self.checkbutton_whole_folder = ttk.Checkbutton(self.output_tab, var=self.run_folder)
        self.checkbutton_whole_folder.grid(column=1, row=11, sticky=W + E, padx=self.padx)
        self.checkbutton_whole_folder.config(state='disabled')

        checkbutton_label = ttk.Label(self.output_tab,
                                      text="Check to generate/run iterations for each file in the directory",
                                      font=self.labelFont)
        checkbutton_label.grid(column=0, row=11, sticky=W)
#___Analysis Tab____________________________________________________________________________________________
    def build_analysis_tab(self):
    
        def select_analysis_folder():
            os.chdir(pathlib.Path.cwd().parent)  # change the working directory from gui to nano-indent

            self.folder_name = filedialog.askdirectory(initialdir=pathlib.Path.cwd(), title="Choose a folder")
            if not self.folder_name:
                self.analysis_dir.set('Please choose a directory')
            else:
                # folder_name = os.path.join(folder_name,'feff')
                self.analysis_dir.set(self.folder_name)

            os.chdir(pathlib.Path.cwd().joinpath('gui'))

        def calculate_and_plot():
            #self.background_types = ['Shirley-Sherwood', 'Slope', 'Exponential', 'Baseline', 'Polynomial 1', 'Polynomial 2', 'Polynomial 3', '3-Param Tougaard', '2-Param Tougaard']

            params = {
                'base': pathlib.Path.cwd().parent,
                'file': self.csv_generate_from,
                'fileName' : self.csv_file.get(),
                'peaks' : self.peak_types,
                'bkgns' : self.background_types,
                'data obj' : self.data_obj


            }

            self.params,self.errors, self.errors_bkgns, self.BTou2, self.BTou3, self.peak_areas, self.FWHM_values, self.peak_y_vals, self.totalFit, self.background_fit, self.residual_fit, self.y_raw, self.upper_error_area, self.lower_error_area = self.analysis_obj.initial_parameters(self.analysis_dir,params,self.scale_var,self.data_KE, self.data_XES, peakType='Voigt',title='Fit')

            #get_params_for_export()
            self.numPeaks = int(self.number_of_peaks.get()) #This only works if you select number of peaks in fitting paramters tab,not if you just want to look at previous fit

            i=int(self.peak_number.get()-1)
            peakType = self.params[i][-1] #Last element in array is the curve fit type

            if i > self.numPeaks:
                print("Number of peaks limit reached")
                pass

            if(peakType.lower() == "lorentzian"):
                BE_text = str(round(self.params[i][0],2)) + "  " + "+/-" + "  " + str(self.errors[i][0])
                self.BE_text = StringVar(value=BE_text)

                lorentzian_text = str(round(self.params[i][1],3)) + "  " + "+/-" + "  " + str(self.errors[i][1])
                self.lorentzian_text = StringVar(value=lorentzian_text)

                amp_text = str(round(self.params[i][2],2)) + "  " + "+/-" + "  " + str(self.errors[i][2])
                self.amp_text = StringVar(value=amp_text)

                area_text = str(round(self.peak_areas[i],2)) + "  " + "+/-" + "  " + str(round(self.upper_error_area[i],3))
                self.area_text = StringVar(value=area_text)

                FWHM_text = str(round(self.FWHM_values[i],2))
                self.FWHM_text = StringVar(value=FWHM_text)

                self.peak_energy_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.BE_text, font=self.entryFont)
                self.peak_energy_entries[i].grid(column=1, row=2, sticky=(W, E))

                self.sigma_entries[i] = ttk.Label(self.analysis_tab, textvariable="0", font=self.entryFont)
                self.sigma_entries[i].grid(column=1, row=3, sticky=(W, E))

                self.lorentzian_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.lorentzian_text, font=self.entryFont)
                self.lorentzian_entries[i].grid(column=1, row=4, sticky=(W, E))

                self.amp_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.amp_text, font=self.entryFont)
                self.amp_entries[i].grid(column=1, row=5, sticky=(W, E))

                self.area_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.area_text, font=self.entryFont)
                self.area_entries[i].grid(column=1, row=6, sticky=(W, E))

                self.FWHM_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.FWHM_text, font=self.entryFont)
                self.FWHM_entries[i].grid(column=1, row=7, sticky=(W, E))


            elif(peakType.lower() == "gaussian"):
                BE_text = str(round(self.params[i][0],2)) + "  " + "+/-" + "  " + str(self.errors[i][0])
                self.BE_text = StringVar(value=BE_text)

                sigma_text = str(round(self.params[i][1],3)) + "  " + "+/-" + "  " + str(self.errors[i][1])
                self.sigma_text = StringVar(value=sigma_text)

                amp_text = str(round(self.params[i][2],2)) + "  " + "+/-" + "  " + str(self.errors[i][2])
                self.amp_text = StringVar(value=amp_text)

                area_text = str(round(self.peak_areas[i],2)) + "  " + "+/-" + "  " + str(round(self.upper_error_area[i],3))
                self.area_text = StringVar(value=area_text)

                FWHM_text = str(round(self.FWHM_values[i],2))
                self.FWHM_text = StringVar(value=FWHM_text)

                self.peak_energy_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.BE_text, font=self.entryFont)
                self.peak_energy_entries[i].grid(column=1, row=2, sticky=(W, E))

                self.sigma_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.sigma_text, font=self.entryFont)
                self.sigma_entries[i].grid(column=1, row=3, sticky=(W, E))

                self.lorentzian_entries[i] = ttk.Label(self.analysis_tab, textvariable="0", font=self.entryFont)
                self.lorentzian_entries[i].grid(column=1, row=4, sticky=(W, E))

                self.amp_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.amp_text, font=self.entryFont)
                self.amp_entries[i].grid(column=1, row=5, sticky=(W, E))

                self.area_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.area_text, font=self.entryFont)
                self.area_entries[i].grid(column=1, row=6, sticky=(W, E))

                self.FWHM_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.FWHM_text, font=self.entryFont)
                self.FWHM_entries[i].grid(column=1, row=7, sticky=(W, E))

            else: #Voigt, Double-Lorentzian, Doniach-Sunjic --> Asymmetry values not being shown right now


                BE_text = str(round(self.params[i][0],2)) + "  " + "+/-" + "  " + str(self.errors[i][0])
                self.BE_text = StringVar(value=BE_text)

                sigma_text = str(round(self.params[i][1],3)) + "  " + "+/-" + "  " + str(self.errors[i][1])
                self.sigma_text = StringVar(value=sigma_text)

                lorentzian_text = str(round(self.params[i][2],3)) + "  " + "+/-" + "  " + str(self.errors[i][2])
                self.lorentzian_text = StringVar(value=lorentzian_text)

                amp_text = str(round(self.params[i][3],2)) + "  " + "+/-" + "  " + str(self.errors[i][3])
                self.amp_text = StringVar(value=amp_text)

                area_text = str(round(self.peak_areas[i],2)) + "  " + "+/-" + "  " + str(round(self.upper_error_area[i],3))
                self.area_text = StringVar(value=area_text)

                FWHM_text = str(round(self.FWHM_values[i],2))
                self.FWHM_text = StringVar(value=FWHM_text)

                self.peak_energy_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.BE_text, font=self.entryFont)
                self.peak_energy_entries[i].grid(column=1, row=2, sticky=(W, E))

                self.sigma_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.sigma_text, font=self.entryFont)
                self.sigma_entries[i].grid(column=1, row=3, sticky=(W, E))

                self.lorentzian_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.lorentzian_text, font=self.entryFont)
                self.lorentzian_entries[i].grid(column=1, row=4, sticky=(W, E))

                self.amp_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.amp_text, font=self.entryFont)
                self.amp_entries[i].grid(column=1, row=5, sticky=(W, E))

                self.area_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.area_text, font=self.entryFont)
                self.area_entries[i].grid(column=1, row=6, sticky=(W, E))

                self.FWHM_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.FWHM_text, font=self.entryFont)
                self.FWHM_entries[i].grid(column=1, row=7, sticky=(W, E))

            #updateExportDataRows(self)

        # def get_params_for_export():
        #     self.analysis_export = self.analysis_obj.get_params(self)
        #     return self.analysis_export

        # Labels
        self.num_peak_labels =[0]*10
        self.sigma_labels = [0]*10
        self.lorentzian_labels = [0]*10
        self.amp_labels = [0]*10
        self.area_labels = [0]*10
        self.FWHM_labels = [0]*10
        # Entries
        self.peak_energy_entries = [0]*10
        self.sigma_entries = [0]*10
        self.lorentzian_entries = [0]*10
        self.amp_entries = [0]*10
        self.area_entries = [0]*10
        self.FWHM_entries = [0]*10
        # Export Data Values
        self.peak_energy_export = [0]*10
        self.sigma_export = [0]*10
        self.lorentzian_export = [0]*10
        self.amp_export = [0]*10
        self.area_export = [0]*10
        self.FWHM_export = [0]*10



        def export_indi_peak_fit():
            x = []
            y = []
            total_fit = []

            x.append(self.data_obj.get_x(self.data_KE, self.data_XES))
            y.append(self.data_obj.get_y(self.scale_var))

            XPS_Oasis_out_name = self.csv_generate_from.stem + "_fit_Aanalyzer" + '.fil'
            os.chdir(pathlib.Path.cwd().parent)
            XPS_Oasis_Fit_folder_path = pathlib.Path(self.folder_name)
            XPS_Oasis_Fit_file = XPS_Oasis_Fit_folder_path.joinpath(XPS_Oasis_out_name)
            os.chdir(pathlib.Path.cwd().joinpath('gui'))
            XPS_Oasis_Fit_file = open(XPS_Oasis_Fit_file, "w")

            XPS_Oasis_Fit_file.write(str("spectraIsIncludedInFil"))
            XPS_Oasis_Fit_file.write(str(" \n"))
            #BE condition:
            XPS_Oasis_Fit_file.write(str("drawBindingEnergy"))
            XPS_Oasis_Fit_file.write(str(" \n"))
            XPS_Oasis_Fit_file.write(str("recoverParameters"))
            XPS_Oasis_Fit_file.write(str(" \n"))
            XPS_Oasis_Fit_file.write(str("tolerance       0.01")) #Idk if we need this parameter in XPS Neo
            XPS_Oasis_Fit_file.write(str(" \n"))
            XPS_Oasis_Fit_file.write(str("instrumentWidth     0.0001")) #Idk if we need this parameter in XPS Neo
            XPS_Oasis_Fit_file.write(str(" \n"))
            XPS_Oasis_Fit_file.write(str("maxNumIterations  4")) # Set to default Aanalyzer value
            XPS_Oasis_Fit_file.write(str(" \n"))
            XPS_Oasis_Fit_file.write(str("iterationsIntegralBkgn  6")) # Set to default Aanalyzer value
            XPS_Oasis_Fit_file.write(str(" \n"))
            XPS_Oasis_Fit_file.write(str("backgroundType"))
            XPS_Oasis_Fit_file.write(str(" \n"))

            #Finding and writing out all the background parameters to fit file
            self.backgrounds_output = []
            self.background_output_values = []
            for i in range(len(self.params)):
                for j in range(len(self.params[i])):
                    if self.params[i][j] == 'Baseline':

                        XPS_Oasis_Fit_file.write(str(" Poly0 poly0Bkgn   "))
                        XPS_Oasis_Fit_file.write(str(self.params[i][j-1]))

                        XPS_Oasis_Fit_file.write(str("  free"))
                        XPS_Oasis_Fit_file.write(str(" \n"))


                    elif self.params[i][j] == 'Polynomial 1':

                        XPS_Oasis_Fit_file.write(str(" Poly1 poly1Bkgn   "))
                        XPS_Oasis_Fit_file.write(str(self.params[i][j-1]))

                        XPS_Oasis_Fit_file.write(str("  free"))
                        XPS_Oasis_Fit_file.write(str(" \n"))

                    elif self.params[i][j] == 'Polynomial 2':

                        XPS_Oasis_Fit_file.write(str(" Poly2 poly2Bkgn   "))
                        XPS_Oasis_Fit_file.write(str(self.params[i][j-1]))

                        XPS_Oasis_Fit_file.write(str("  free"))
                        XPS_Oasis_Fit_file.write(str(" \n"))

                    elif self.params[i][j] == 'Polynomial 3':

                        XPS_Oasis_Fit_file.write(str(" Poly3 poly3Bkgn   "))
                        XPS_Oasis_Fit_file.write(str(self.params[i][j-1]))

                        XPS_Oasis_Fit_file.write(str("  free"))
                        XPS_Oasis_Fit_file.write(str(" \n"))

                    elif self.params[i][j] == 'Exponential': #Exponential is done in combo with tougaard CTou2PActive/CTou3PActive line IDK if it will work if not written like that

                        XPS_Oasis_Fit_file.write(str(" Exp exponent   "))
                        XPS_Oasis_Fit_file.write(str(self.params[i][j-1]))
                        XPS_Oasis_Fit_file.write(str("  free"))
                        XPS_Oasis_Fit_file.write(str(" \n"))
                    
                    elif self.params[i][j] == 'SVSC': #changed from Shirley-Sherwood
                        self.backgrounds_output.append(self.params[i][j])
                        self.background_output_values.append(self.params[i][j-1])
                        '''
                        XPS_Oasis_Fit_file.write(str(" Int integralBkgn ")) #IDK why but Aanalyzer has it that this is called something else when Slope AND Shirley are selected
                        XPS_Oasis_Fit_file.write(str(self.params[i][j-1])) #IntSlope integralSlopeBkgn
                        XPS_Oasis_Fit_file.write(str("  free"))
                        XPS_Oasis_Fit_file.write(str(" \n"))
                    '''

                    elif self.params[i][j] == 'Shirley':


                        XPS_Oasis_Fit_file.write(str(" Int integralBkgn ")) #IDK why but Aanalyzer has it that this is called something else when Slope AND Shirley are selected
                        XPS_Oasis_Fit_file.write(str(self.params[i][j-1])) #IntSlope integralSlopeBkgn
                        XPS_Oasis_Fit_file.write(str("  free"))
                        XPS_Oasis_Fit_file.write(str(" \n"))

                    elif self.params[i][j] == 'Linear':


                        XPS_Oasis_Fit_file.write(str(" Int integralBkgn "))
                        XPS_Oasis_Fit_file.write(str(self.params[i][j-1]))
                        XPS_Oasis_Fit_file.write(str("  free"))
                        XPS_Oasis_Fit_file.write(str(" \n"))


                    elif self.params[i][j] == '2-Param Tougaard':


                        XPS_Oasis_Fit_file.write(str(" useTou2PActive BTou2PActive      "))
                        XPS_Oasis_Fit_file.write(str(round(self.BTou2, 2))) #We need to have this parameter printed out #B2
                        XPS_Oasis_Fit_file.write(str("  fix")) #This is currently the only parameter being held fixed
                        XPS_Oasis_Fit_file.write(str(" \n"))
                        XPS_Oasis_Fit_file.write(str(" CTou2PActive      "))
                        XPS_Oasis_Fit_file.write(str(1643))
                        XPS_Oasis_Fit_file.write(str(" \n"))

                    elif self.params[i][j] == '3-Param Tougaard':


                        XPS_Oasis_Fit_file.write(str(" useTou3PActive BTou3PActive      "))
                        XPS_Oasis_Fit_file.write(str(round(self.BTou3, 2))) #We need to have this parameter printed out #B3
                        XPS_Oasis_Fit_file.write(str("  free"))
                        XPS_Oasis_Fit_file.write(str(" \n"))
                        XPS_Oasis_Fit_file.write(str(" CTou3PActive      "))
                        XPS_Oasis_Fit_file.write(str(1000))
                        XPS_Oasis_Fit_file.write(str("    "))
                        XPS_Oasis_Fit_file.write(str(" DTou3PActive      "))
                        XPS_Oasis_Fit_file.write(str(13300))
                        XPS_Oasis_Fit_file.write(str("    "))
                        XPS_Oasis_Fit_file.write(str(" \n"))

            XPS_Oasis_Fit_file.write(str("  doNotAllowForNegativeSVSC"))
            XPS_Oasis_Fit_file.write(str(" \n"))
            XPS_Oasis_Fit_file.write(str("nextData 0"))
            XPS_Oasis_Fit_file.write(str(" \n"))
            XPS_Oasis_Fit_file.write(str("fileName"))
            XPS_Oasis_Fit_file.write(str(" "))
            XPS_Oasis_Fit_file.write(str(self.csv_generate_from.stem))
            XPS_Oasis_Fit_file.write(str(".txt"))
            XPS_Oasis_Fit_file.write(str(" \n"))
            XPS_Oasis_Fit_file.write(str(" photonEnergy     1486.7 dataSpectrumStarts")) #Even if this is different we already give the enegry in BE so this doesnt matter. It does if we give it in KE...fix
            XPS_Oasis_Fit_file.write(str(" \n"))







            #Changed output format to be similiar to Aanalyzer. Only have to output one file to enter into Igor pro to make graphs.
            peak_num = len(self.peak_y_vals)
            out_name = self.csv_generate_from.stem + "_Fit_Arrays" + '.txt'
            os.chdir(pathlib.Path.cwd().parent)
            peakFit_folder_path = pathlib.Path(self.folder_name)
            peakFit_file = peakFit_folder_path.joinpath(out_name)
            os.chdir(pathlib.Path.cwd().joinpath('gui'))
            peakFit_file = open(peakFit_file, "w")
            peakFit_file.write(str("x"))
            peakFit_file.write(str("  "))
            peakFit_file.write(str("yExp"))
            peakFit_file.write(str("  "))
            peakFit_file.write(str("yCal"))
            peakFit_file.write(str("  "))
            peakFit_file.write(str("Bkgn"))
            peakFit_file.write(str("  "))
            peakFit_file.write(str("Res"))
            for i in range(peak_num):
                peakFit_file.write(str("  "))
                peakFit_file.write(str("peak")+str(i+1))
            peakFit_file.write(str(" \n"))
            for j in range(len(x[0])):
                peakFit_file.write(str(x[0][j]))
                peakFit_file.write(str("  "))
                peakFit_file.write(str(self.y_raw[j]))
                peakFit_file.write(str("  "))
                peakFit_file.write(str(self.totalFit[j]))
                peakFit_file.write(str("  "))
                peakFit_file.write(str(self.background_fit[j]))
                peakFit_file.write(str("  "))
                peakFit_file.write(str(self.totalFit[j]-self.y_raw[j]))
                peakFit_file.write(str("  "))
                for i in range(peak_num):
                    peakFit_file.write(str(self.peak_y_vals[i][j]))
                    peakFit_file.write(str("  "))
                peakFit_file.write(str(" \n"))
      
            
            peakFit_file.close()













            '''
            for i in range(len(self.peak_y_vals)):

                out_name = self.csv_generate_from.stem + "_" + "Peak" + "_" + str(i+1) + "_" + "fit" + '.txt'
                total_out_name = self.csv_generate_from.stem + "_" + "Total" + "_" + "fit" + '.txt'
                background_out_name = self.csv_generate_from.stem + "_" + "Background" + "_" + "fit" + '.txt'
                residual_out_name = self.csv_generate_from.stem + "_" + "Residual" + "_" + "fit" + '.txt'
                os.chdir(pathlib.Path.cwd().parent)


                peakFit_folder_path = pathlib.Path(self.folder_name)
                peakFit_file = peakFit_folder_path.joinpath(out_name)
                totalPeakFit_folder_path = pathlib.Path(self.folder_name)
                totalPeakFit_file = totalPeakFit_folder_path.joinpath(total_out_name)
                backgroundFit_folder_path = pathlib.Path(self.folder_name)
                backgroundFit_file = backgroundFit_folder_path.joinpath(background_out_name)
                residualFit_folder_path = pathlib.Path(self.folder_name)
                residualFit_file = residualFit_folder_path.joinpath(residual_out_name)

                os.chdir(pathlib.Path.cwd().joinpath('gui'))



                peakFit_file = open(peakFit_file, "w")
                totalPeakFit_file = open(totalPeakFit_file, "w")
                backgroundFit_file = open(backgroundFit_file, "w")
                residualFit_file = open(residualFit_file, "w")


                peakFit_file.write(str("x"))
                peakFit_file.write(str("  "))
                peakFit_file.write(str("yCal"))
                peakFit_file.write(str("  "))
                peakFit_file.write(str("Bkgn"))
                peakFit_file.write(str(" \n"))
                for j in range(len(x[0])):

                    peakFit_file.write(str(x[0][j]))
                    peakFit_file.write(str("  "))
                    peakFit_file.write(str(self.peak_y_vals[i][j])) #Individual peak fit components
                    peakFit_file.write(str(" \n"))
                    totalPeakFit_file.write(str(x[0][j]))
                    totalPeakFit_file.write(str("  "))
                    totalPeakFit_file.write(str(self.totalFit[j]))
                    totalPeakFit_file.write(str(" \n"))
                    backgroundFit_file.write(str(x[0][j]))
                    backgroundFit_file.write(str("  "))
                    backgroundFit_file.write(str(self.background_fit[j]))
                    backgroundFit_file.write(str(" \n"))
                    residualFit_file.write(str(x[0][j]))
                    residualFit_file.write(str("  "))
                    residualFit_file.write(str(self.residual_fit[j]))
                    residualFit_file.write(str(" \n"))


                peakFit_file.close()
                totalPeakFit_file.close()
                backgroundFit_file.close()

                print("Peak ", i+1, " information saved to output file location")
            '''
            photonEnergy = 1486.6


            for j in range(len(x[0])):
                if x[0][0] < x[0][-1]: #KE condition
                    x_value = x[0][j]
                else:
                    x_value = photonEnergy - x[0][j]
                XPS_Oasis_Fit_file.write(str(x_value))
                XPS_Oasis_Fit_file.write(str("  "))
                XPS_Oasis_Fit_file.write(str(y[0][j]))
                XPS_Oasis_Fit_file.write(str(" \n"))

            XPS_Oasis_Fit_file.write(str("dataSpectrumEnds"))
            XPS_Oasis_Fit_file.write(str(" \n"))
            XPS_Oasis_Fit_file.write(str("comment")) #Allow for users to input comments?
            XPS_Oasis_Fit_file.write(str(" \n"))
            XPS_Oasis_Fit_file.write(str("skipLines"))
            XPS_Oasis_Fit_file.write(str(" "))
            XPS_Oasis_Fit_file.write(str(self.skipLn.get()))
            XPS_Oasis_Fit_file.write(str(" \n"))
            XPS_Oasis_Fit_file.write(str("xChannel 0 yChannel 1"))
            XPS_Oasis_Fit_file.write(str(" \n"))
            XPS_Oasis_Fit_file.write(str(" shift         0  offset         0  gain         1  externalX         0"))
            XPS_Oasis_Fit_file.write(str(" \n"))
            XPS_Oasis_Fit_file.write(str("Parameters"))
            XPS_Oasis_Fit_file.write(str(" \n"))
            #XPS_Oasis_Fit_file.write(str("recoverParameters"))
            #XPS_Oasis_Fit_file.write(str(" \n"))
            #XPS_Oasis_Fit_file.write(str("leftLimit"))
            #XPS_Oasis_Fit_file.write(str("    "))
            #XPS_Oasis_Fit_file.write(str(x[0][0]))
            #XPS_Oasis_Fit_file.write(str("   "))
            #XPS_Oasis_Fit_file.write(str("rightLimit"))
            #XPS_Oasis_Fit_file.write(str("    "))
            #XPS_Oasis_Fit_file.write(str(x[0][-1]))
            #XPS_Oasis_Fit_file.write(str(" \n"))
            

            #Getting fit components of each peak
            ck_count = 0
    
            for k in range(len(self.peak_y_vals)):
                XPS_Oasis_Fit_file.write(str("nextPeak "))
                XPS_Oasis_Fit_file.write(str(k + ck_count))
                XPS_Oasis_Fit_file.write(str(" \n"))


               


                peakType = self.params[k][-1]

                photonEnergy = 1486.6 #Want as user input later
                if x[0][0] < x[0][-1]: #KE condition
                    KE = self.params[k][0]
                else:
                    KE = photonEnergy - self.params[k][0] #Values are read in BE in Aanalyzer


               
                if(peakType.lower() == "lorentzian"):
                    if len(self.params[k]) >= 5:
                        XPS_Oasis_Fit_file.write(str("doublet"))
                        XPS_Oasis_Fit_file.write(str(" \n"))
                        XPS_Oasis_Fit_file.write(str("soSplitting  "))
                        XPS_Oasis_Fit_file.write(str(self.params[k][3])) #This needs to be negative in Aanalyzer
                        XPS_Oasis_Fit_file.write(str(" free"))
                        XPS_Oasis_Fit_file.write(str(" \n"))
                        XPS_Oasis_Fit_file.write(str("branchingRatio      "))
                        XPS_Oasis_Fit_file.write(str("0.75")) #WE DO NOT HAVE A WAY OF OUTPUTTING THIS RIGHT NOW --> FIX THIS
                        XPS_Oasis_Fit_file.write(str(" fix"))
                        XPS_Oasis_Fit_file.write(str(" \n"))

                    else:
                        XPS_Oasis_Fit_file.write(str("singlet"))
                        XPS_Oasis_Fit_file.write(str(" \n"))

                    XPS_Oasis_Fit_file.write(str("gauss   1 limited"))
                    XPS_Oasis_Fit_file.write(str(" \n"))
                    XPS_Oasis_Fit_file.write(str("lorentzian      "))
                    XPS_Oasis_Fit_file.write(str(self.params[k][1]))
                    XPS_Oasis_Fit_file.write(str(" free"))
                    XPS_Oasis_Fit_file.write(str(" \n"))
                    XPS_Oasis_Fit_file.write(str("energy   "))
                    XPS_Oasis_Fit_file.write(str(KE))
                    XPS_Oasis_Fit_file.write(str(" free"))
                    XPS_Oasis_Fit_file.write(str(" \n"))
                    XPS_Oasis_Fit_file.write(str("asymmetry         1  limitUp    50 limitDown     1"))
                    XPS_Oasis_Fit_file.write(str(" \n"))
                    XPS_Oasis_Fit_file.write(str("asymmetryDoniach       0.1  limitUp     1 limitDown     0"))
                    XPS_Oasis_Fit_file.write(str(" \n"))
                    XPS_Oasis_Fit_file.write(str("area   "))
                    XPS_Oasis_Fit_file.write(str(self.peak_areas[k]))
                    XPS_Oasis_Fit_file.write(str(" free"))
                    XPS_Oasis_Fit_file.write(str(" \n"))
                    XPS_Oasis_Fit_file.write(str("peakType "))
                    XPS_Oasis_Fit_file.write(str(peakType))
                    XPS_Oasis_Fit_file.write(str(" \n"))


                    for i in range(len(self.params)):
                        for j in range(len(self.params[i])):
                            if self.params[i][j] == 'SVSC': #changed from Shirley-Sherwood
                                XPS_Oasis_Fit_file.write(str("peakShirleyBackground peakShirleyBackgroundValue      ")) 
                                XPS_Oasis_Fit_file.write(str(self.background_output_values[k])) #IntSlope integralSlopeBkgn
                                XPS_Oasis_Fit_file.write(str("  fix"))
                                XPS_Oasis_Fit_file.write(str(" \n"))
                                break

                    XPS_Oasis_Fit_file.write(str("peakColor clOlive")) #Change this based on k
                    XPS_Oasis_Fit_file.write(str(" \n"))

                elif(peakType.lower() == "gaussian"):
                    if len(self.params[k]) >= 5:
                        XPS_Oasis_Fit_file.write(str("doublet"))
                        XPS_Oasis_Fit_file.write(str(" \n"))
                        XPS_Oasis_Fit_file.write(str("soSplitting  "))
                        XPS_Oasis_Fit_file.write(str(self.params[k][3]))
                        XPS_Oasis_Fit_file.write(str(" free"))
                        XPS_Oasis_Fit_file.write(str(" \n"))
                        XPS_Oasis_Fit_file.write(str("branchingRatio      "))
                        XPS_Oasis_Fit_file.write(str("0.75")) #WE DO NOT HAVE A WAY OF OUTPUTTING THIS RIGHT NOW --> FIX THIS
                        XPS_Oasis_Fit_file.write(str(" fix"))
                        XPS_Oasis_Fit_file.write(str(" \n"))
                    else:
                        XPS_Oasis_Fit_file.write(str("singlet"))
                        XPS_Oasis_Fit_file.write(str(" \n"))
                    XPS_Oasis_Fit_file.write(str("gauss   "))
                    XPS_Oasis_Fit_file.write(str(self.params[k][1]))
                    XPS_Oasis_Fit_file.write(str(" free"))
                    XPS_Oasis_Fit_file.write(str(" \n"))
                    XPS_Oasis_Fit_file.write(str("lorentzian      0.27 limited"))
                    XPS_Oasis_Fit_file.write(str(" \n"))
                    XPS_Oasis_Fit_file.write(str("energy   "))
                    XPS_Oasis_Fit_file.write(str(KE))
                    XPS_Oasis_Fit_file.write(str(" free"))
                    XPS_Oasis_Fit_file.write(str(" \n"))
                    XPS_Oasis_Fit_file.write(str("asymmetry         1  limitUp    50 limitDown     1"))
                    XPS_Oasis_Fit_file.write(str(" \n"))
                    XPS_Oasis_Fit_file.write(str("asymmetryDoniach       0.1  limitUp     1 limitDown     0"))
                    XPS_Oasis_Fit_file.write(str(" \n"))
                    XPS_Oasis_Fit_file.write(str("area   "))
                    XPS_Oasis_Fit_file.write(str(self.peak_areas[k]))
                    XPS_Oasis_Fit_file.write(str(" free"))
                    XPS_Oasis_Fit_file.write(str(" \n"))
                    XPS_Oasis_Fit_file.write(str("peakType "))
                    XPS_Oasis_Fit_file.write(str(peakType))
                    XPS_Oasis_Fit_file.write(str(" \n"))


                    for i in range(len(self.params)):
                        for j in range(len(self.params[i])):
                            if self.params[i][j] == 'SVSC': #changed from Shirley-Sherwood
                                XPS_Oasis_Fit_file.write(str("peakShirleyBackground peakShirleyBackgroundValue      ")) 
                                XPS_Oasis_Fit_file.write(str(self.background_output_values[k])) #IntSlope integralSlopeBkgn
                                XPS_Oasis_Fit_file.write(str("  fix"))
                                XPS_Oasis_Fit_file.write(str(" \n"))
                                break

                    XPS_Oasis_Fit_file.write(str("peakColor clOlive")) #Change this based on k
                    XPS_Oasis_Fit_file.write(str(" \n"))

                elif(peakType.lower() == "voigt"):
                    if len(self.params[k]) >= 6:
                        XPS_Oasis_Fit_file.write(str("doublet"))
                        XPS_Oasis_Fit_file.write(str(" \n"))
                        XPS_Oasis_Fit_file.write(str("soSplitting  "))
                        XPS_Oasis_Fit_file.write(str(-self.params[k][5]))
                        XPS_Oasis_Fit_file.write(str(" free"))
                        XPS_Oasis_Fit_file.write(str(" \n"))
                        XPS_Oasis_Fit_file.write(str("branchingRatio      "))
                        XPS_Oasis_Fit_file.write(str(self.params[k][4])) #WE DO NOT HAVE A WAY OF OUTPUTTING THIS RIGHT NOW --> FIX THIS
                        XPS_Oasis_Fit_file.write(str(" fix"))
                        XPS_Oasis_Fit_file.write(str(" \n"))
                    else:
                        XPS_Oasis_Fit_file.write(str("singlet"))
                        XPS_Oasis_Fit_file.write(str(" \n"))
                    XPS_Oasis_Fit_file.write(str("gauss   "))
                    XPS_Oasis_Fit_file.write(str(self.params[k][1]))
                    XPS_Oasis_Fit_file.write(str(" free"))
                    XPS_Oasis_Fit_file.write(str(" \n"))
                    XPS_Oasis_Fit_file.write(str("lorentzian      "))
                    XPS_Oasis_Fit_file.write(str(self.params[k][2]))
                    XPS_Oasis_Fit_file.write(str(" free"))
                    XPS_Oasis_Fit_file.write(str(" \n"))
                    XPS_Oasis_Fit_file.write(str("energy   "))
                    XPS_Oasis_Fit_file.write(str(KE))
                    XPS_Oasis_Fit_file.write(str(" free"))
                    XPS_Oasis_Fit_file.write(str(" \n"))
                    XPS_Oasis_Fit_file.write(str("asymmetry         1  limitUp    50 limitDown     1"))
                    XPS_Oasis_Fit_file.write(str(" \n"))
                    XPS_Oasis_Fit_file.write(str("asymmetryDoniach       0.1  limitUp     1 limitDown     0"))
                    XPS_Oasis_Fit_file.write(str(" \n"))
                    XPS_Oasis_Fit_file.write(str("area   "))
                    XPS_Oasis_Fit_file.write(str(self.peak_areas[k]))
                    XPS_Oasis_Fit_file.write(str(" free"))
                    XPS_Oasis_Fit_file.write(str(" \n"))
                    XPS_Oasis_Fit_file.write(str("peakType "))
                    XPS_Oasis_Fit_file.write(str(peakType))
                    XPS_Oasis_Fit_file.write(str(" \n"))



                    for i in range(len(self.params)):
                        for j in range(len(self.params[i])):
                            if self.params[i][j] == 'SVSC': #changed from Shirley-Sherwood
                                XPS_Oasis_Fit_file.write(str("peakShirleyBackground peakShirleyBackgroundValue      ")) 
                                XPS_Oasis_Fit_file.write(str(self.background_output_values[k])) #IntSlope integralSlopeBkgn
                                XPS_Oasis_Fit_file.write(str("  fix"))
                                XPS_Oasis_Fit_file.write(str(" \n"))
                                break

                    XPS_Oasis_Fit_file.write(str("peakColor clOlive")) #Change this based on k
                    XPS_Oasis_Fit_file.write(str(" \n"))

                elif(peakType.lower() == "double lorentzian"):
                    if len(self.params[k]) >= 6:

                        if len(self.params[k]) >= 10: #If Coster-Kronig
                            XPS_Oasis_Fit_file.write(str("singlet"))
                            XPS_Oasis_Fit_file.write(str(" \n"))
                            XPS_Oasis_Fit_file.write(str("gauss   "))
                            XPS_Oasis_Fit_file.write(str(self.params[k][1]))
                            XPS_Oasis_Fit_file.write(str(" free"))
                            XPS_Oasis_Fit_file.write(str(" \n"))
                            XPS_Oasis_Fit_file.write(str("lorentzian      "))
                            XPS_Oasis_Fit_file.write(str(self.params[k][2]))
                            XPS_Oasis_Fit_file.write(str(" free"))
                            XPS_Oasis_Fit_file.write(str(" \n"))
                            XPS_Oasis_Fit_file.write(str("energy   "))
                            XPS_Oasis_Fit_file.write(str(KE))
                            XPS_Oasis_Fit_file.write(str(" free"))
                            XPS_Oasis_Fit_file.write(str(" \n"))
                            XPS_Oasis_Fit_file.write(str("asymmetry         "))
                            XPS_Oasis_Fit_file.write(str(self.params[k][4]))
                            XPS_Oasis_Fit_file.write(str("  limitUp    50 limitDown     1"))
                            XPS_Oasis_Fit_file.write(str(" \n"))
                            XPS_Oasis_Fit_file.write(str("asymmetryDoniach       0.1  limitUp     1 limitDown     0"))
                            XPS_Oasis_Fit_file.write(str(" \n"))
                            XPS_Oasis_Fit_file.write(str("area   "))
                            XPS_Oasis_Fit_file.write(str((self.peak_areas[k])/(1+self.params[k][5]))) #this is the total area, we need to limit this to correspond to two peaks
                            XPS_Oasis_Fit_file.write(str(" free"))
                            XPS_Oasis_Fit_file.write(str(" \n"))
                            XPS_Oasis_Fit_file.write(str("peakType "))
                            XPS_Oasis_Fit_file.write(str("DoubleLorentzian"))
                            XPS_Oasis_Fit_file.write(str(" \n"))



                            for i in range(len(self.params)):
                                for j in range(len(self.params[i])):
                                    if self.params[i][j] == 'SVSC': #changed from Shirley-Sherwood
                                        XPS_Oasis_Fit_file.write(str("peakShirleyBackground peakShirleyBackgroundValue      ")) 
                                        XPS_Oasis_Fit_file.write(str(self.background_output_values[k])) #IntSlope integralSlopeBkgn
                                        XPS_Oasis_Fit_file.write(str("  fix"))
                                        XPS_Oasis_Fit_file.write(str(" \n"))
                                        break

                            
                            XPS_Oasis_Fit_file.write(str("peakColor clOlive")) #Change this based on k
                            XPS_Oasis_Fit_file.write(str(" \n"))


                            XPS_Oasis_Fit_file.write(str("nextPeak ")) #How to account for more peaks after this??? add count += 1 to it?
                            ck_count += 1
                            XPS_Oasis_Fit_file.write(str(k+ck_count))
                            XPS_Oasis_Fit_file.write(str(" \n"))


                            XPS_Oasis_Fit_file.write(str("singlet"))
                            XPS_Oasis_Fit_file.write(str(" \n"))
                            XPS_Oasis_Fit_file.write(str("gauss   "))
                            #gauss = self.params[k][1] * 2 * np.sqrt(2*np.log(2.0)) # dont think we need this anymore? CHECK
                            XPS_Oasis_Fit_file.write(str(self.params[k][1]))
                            XPS_Oasis_Fit_file.write(str(" correlation g")) #Need to get peak not just 0 CHECK
                            XPS_Oasis_Fit_file.write(str(k)) #which peak we correlate to
                            XPS_Oasis_Fit_file.write(str(" \n"))
                            XPS_Oasis_Fit_file.write(str("lorentzian      "))
                            XPS_Oasis_Fit_file.write(str(self.params[k][7]))
                            XPS_Oasis_Fit_file.write(str(" free"))
                            XPS_Oasis_Fit_file.write(str(" \n"))
                            XPS_Oasis_Fit_file.write(str("energy   "))
                            XPS_Oasis_Fit_file.write(str(KE))
                            XPS_Oasis_Fit_file.write(str(" correlation e"))
                            XPS_Oasis_Fit_file.write(str(k)) #which peak we correlate to
                            XPS_Oasis_Fit_file.write(str(-self.params[k][6]))
                            XPS_Oasis_Fit_file.write(str(" \n"))
                            XPS_Oasis_Fit_file.write(str("asymmetry         "))
                            XPS_Oasis_Fit_file.write(str(self.params[k][8]))
                            XPS_Oasis_Fit_file.write(str("  limitUp    50 limitDown     1"))
                            XPS_Oasis_Fit_file.write(str(" \n"))
                            XPS_Oasis_Fit_file.write(str("asymmetryDoniach       0.1  limitUp     1 limitDown     0"))
                            XPS_Oasis_Fit_file.write(str(" \n"))
                            XPS_Oasis_Fit_file.write(str("area   "))
                            XPS_Oasis_Fit_file.write(str(((self.peak_areas[k])/(1+self.params[k][5]))*self.params[k][5])) #a1 * Br
                            XPS_Oasis_Fit_file.write(str(" correlation a"))
                            XPS_Oasis_Fit_file.write(str(k)) #which peak we correlate to
                            XPS_Oasis_Fit_file.write(str("*"))
                            XPS_Oasis_Fit_file.write(str(self.params[k][5]))
                            XPS_Oasis_Fit_file.write(str(" \n"))
                            XPS_Oasis_Fit_file.write(str("peakType "))
                            XPS_Oasis_Fit_file.write(str("DoubleLorentzian"))
                            XPS_Oasis_Fit_file.write(str(" \n"))


                            for i in range(len(self.params)):
                                for j in range(len(self.params[i])):
                                    if self.params[i][j] == 'SVSC': #changed from Shirley-Sherwood
                                        XPS_Oasis_Fit_file.write(str("peakShirleyBackground peakShirleyBackgroundValue      ")) 
                                        XPS_Oasis_Fit_file.write(str(self.background_output_values[k])) #IntSlope integralSlopeBkgn
                                        XPS_Oasis_Fit_file.write(str("  fix"))
                                        XPS_Oasis_Fit_file.write(str(" \n"))
                                        break

                            XPS_Oasis_Fit_file.write(str("peakColor clOlive")) #Change this based on k
                            XPS_Oasis_Fit_file.write(str(" \n"))



                        else:    
                            XPS_Oasis_Fit_file.write(str("doublet"))
                            XPS_Oasis_Fit_file.write(str(" \n"))
                            XPS_Oasis_Fit_file.write(str("soSplitting  "))
                            XPS_Oasis_Fit_file.write(str(-self.params[k][6]))
                            XPS_Oasis_Fit_file.write(str(" free"))
                            XPS_Oasis_Fit_file.write(str(" \n"))
                            XPS_Oasis_Fit_file.write(str("branchingRatio      "))
                            XPS_Oasis_Fit_file.write(str(self.params[k][5])) #WE DO NOT HAVE A WAY OF OUTPUTTING THIS RIGHT NOW --> FIX THIS
                            XPS_Oasis_Fit_file.write(str(" fix"))
                            XPS_Oasis_Fit_file.write(str(" \n"))
                            XPS_Oasis_Fit_file.write(str("gauss   "))
                            XPS_Oasis_Fit_file.write(str(self.params[k][1]))
                            XPS_Oasis_Fit_file.write(str(" free"))
                            XPS_Oasis_Fit_file.write(str(" \n"))
                            XPS_Oasis_Fit_file.write(str("lorentzian      "))
                            XPS_Oasis_Fit_file.write(str(self.params[k][2]))
                            XPS_Oasis_Fit_file.write(str(" free"))
                            XPS_Oasis_Fit_file.write(str(" \n"))
                            XPS_Oasis_Fit_file.write(str("energy   "))
                            XPS_Oasis_Fit_file.write(str(KE))
                            XPS_Oasis_Fit_file.write(str(" free"))
                            XPS_Oasis_Fit_file.write(str(" \n"))
                            XPS_Oasis_Fit_file.write(str("asymmetry         "))
                            XPS_Oasis_Fit_file.write(str(self.params[k][4]))
                            XPS_Oasis_Fit_file.write(str("  limitUp    50 limitDown     1"))
                            XPS_Oasis_Fit_file.write(str(" \n"))
                            XPS_Oasis_Fit_file.write(str("asymmetryDoniach       0.1  limitUp     1 limitDown     0"))
                            XPS_Oasis_Fit_file.write(str(" \n"))
                            XPS_Oasis_Fit_file.write(str("area   "))
                            XPS_Oasis_Fit_file.write(str(self.peak_areas[k]))
                            XPS_Oasis_Fit_file.write(str(" free"))
                            XPS_Oasis_Fit_file.write(str(" \n"))
                            XPS_Oasis_Fit_file.write(str("peakType "))
                            XPS_Oasis_Fit_file.write(str("DoubleLorentzian"))
                            XPS_Oasis_Fit_file.write(str(" \n"))


                            for i in range(len(self.params)):
                                for j in range(len(self.params[i])):
                                    if self.params[i][j] == 'SVSC': #changed from Shirley-Sherwood
                                        XPS_Oasis_Fit_file.write(str("peakShirleyBackground peakShirleyBackgroundValue      ")) 
                                        XPS_Oasis_Fit_file.write(str(self.background_output_values[k])) #IntSlope integralSlopeBkgn
                                        XPS_Oasis_Fit_file.write(str("  fix"))
                                        XPS_Oasis_Fit_file.write(str(" \n"))
                                        break


                            XPS_Oasis_Fit_file.write(str("peakColor clOlive")) #Change this based on k
                            XPS_Oasis_Fit_file.write(str(" \n"))
                    else:
                        XPS_Oasis_Fit_file.write(str("singlet"))
                        XPS_Oasis_Fit_file.write(str(" \n"))
                        XPS_Oasis_Fit_file.write(str("gauss   "))
                        XPS_Oasis_Fit_file.write(str(self.params[k][1]))
                        XPS_Oasis_Fit_file.write(str(" free"))
                        XPS_Oasis_Fit_file.write(str(" \n"))
                        XPS_Oasis_Fit_file.write(str("lorentzian      "))
                        XPS_Oasis_Fit_file.write(str(self.params[k][2]))
                        XPS_Oasis_Fit_file.write(str(" free"))
                        XPS_Oasis_Fit_file.write(str(" \n"))
                        XPS_Oasis_Fit_file.write(str("energy   "))
                        XPS_Oasis_Fit_file.write(str(KE))
                        XPS_Oasis_Fit_file.write(str(" free"))
                        XPS_Oasis_Fit_file.write(str(" \n"))
                        XPS_Oasis_Fit_file.write(str("asymmetry         "))
                        XPS_Oasis_Fit_file.write(str(self.params[k][4]))
                        XPS_Oasis_Fit_file.write(str("  limitUp    50 limitDown     1"))
                        XPS_Oasis_Fit_file.write(str(" \n"))
                        XPS_Oasis_Fit_file.write(str("asymmetryDoniach       0.1  limitUp     1 limitDown     0"))
                        XPS_Oasis_Fit_file.write(str(" \n"))
                        XPS_Oasis_Fit_file.write(str("area   "))
                        XPS_Oasis_Fit_file.write(str(self.peak_areas[k]))
                        XPS_Oasis_Fit_file.write(str(" free"))
                        XPS_Oasis_Fit_file.write(str(" \n"))
                        XPS_Oasis_Fit_file.write(str("peakType "))
                        XPS_Oasis_Fit_file.write(str("DoubleLorentzian"))
                        XPS_Oasis_Fit_file.write(str(" \n"))

                        for i in range(len(self.params)):
                            for j in range(len(self.params[i])):
                                if self.params[i][j] == 'SVSC': #changed from Shirley-Sherwood
                                    XPS_Oasis_Fit_file.write(str("peakShirleyBackground peakShirleyBackgroundValue      ")) 
                                    XPS_Oasis_Fit_file.write(str(self.background_output_values[k])) #IntSlope integralSlopeBkgn
                                    XPS_Oasis_Fit_file.write(str("  fix"))
                                    XPS_Oasis_Fit_file.write(str(" \n"))
                                    break

                        XPS_Oasis_Fit_file.write(str("peakColor clOlive")) #Change this based on k
                        XPS_Oasis_Fit_file.write(str(" \n"))


            XPS_Oasis_Fit_file.write(str("  endOfParameters"))
            XPS_Oasis_Fit_file.write(str(" \n"))






            XPS_Oasis_Fit_file.close()
            print("Fit saved to output file location")
            





        #Function for exporting peak parameters as LaTeX table:
        def export_peak_parameters():




            out_name = self.csv_generate_from.stem + "_" + "LaTeX_Table" + '.txt'

            os.chdir(pathlib.Path.cwd().parent)  


            table_folder_path = pathlib.Path(self.folder_name)
            table_file = table_folder_path.joinpath(out_name)
            os.chdir(pathlib.Path.cwd().joinpath('gui'))

            #name_of_file = str("Output_table")

            #completeName = os.path.abspath(name_of_file+".txt")

            table_file = open(table_file, "w")

            #Header and peak rows
            headers = []
            BE_vals = []
            sigma_vals = []
            gamma_vals = []
            amp_vals = []
            area_vals = []
            fwhm_vals = []

            #DL, doublet and CK stuff
            asym_vals = []
            Br_vals = []
            SOS_vals = []
            CK_vals = []

            doublet = False
            asym = False
            CK = False

            data = dict()
            
            headers.append("Peaks") #Error here for some reason. "invalid escape error"
            peak_range = len(self.params) - len(self.errors_bkgns) #To save first press "Plot Best Fit"
            
            for i in range(len(self.peak_y_vals)):
            #for i in range(0, peak_range): #appending the values of each peak to each specified row in the table

                headers.append(str(i+1))
                peakType = self.params[i][-1] #Last element in array is the curve fit type
                is_singlet = self.params[i][-2]
               

                if(peakType.lower() == "lorentzian"):
                    BE_vals.append(str(round(self.params[i][0],2)) + " \\textpm " +  str(self.errors[i][0]))
                    sigma_vals.append(" ")
                    gamma_vals.append(str(round(self.params[i][1],3)) + " \\textpm " +  str(self.errors[i][1]))
                    amp_vals.append(str(round(self.params[i][2],2)) + " \\textpm " +  str(self.errors[i][2]))
                    area_vals.append(str(round(self.peak_areas[i],2)) + " \\textpm " +  str(round(self.upper_error_area[i],3)))
                    fwhm_vals.append(str(round(self.FWHM_values[i],2)))
                elif(peakType.lower() == "gaussian"):
                    BE_vals.append(str(round(self.params[i][0],2)) + " \\textpm " +  str(self.errors[i][0]))
                    sigma_vals.append(str(round(self.params[i][1],3)) + " \\textpm " +  str(self.errors[i][1]))
                    gamma_vals.append(" ")
                    amp_vals.append(str(round(self.params[i][2],2)) + " \\textpm " +  str(self.errors[i][2]))
                    area_vals.append(str(round(self.peak_areas[i],2)) + " \\textpm " +  str(round(self.upper_error_area[i],3)))
                    fwhm_vals.append(str(round(self.FWHM_values[i],2)))
                else: #Does not differentiate between Voigt or Double Lorentzian or Doniach Sunjic
                    BE_vals.append(str(round(self.params[i][0],2)) + " \\textpm " +  str(self.errors[i][0]))
                    sigma_vals.append(str(round(self.params[i][1],3)) + " \\textpm " +  str(self.errors[i][1]))
                    gamma_vals.append(str(round(self.params[i][2],3)) + " \\textpm " +  str(self.errors[i][2]))
                    amp_vals.append(str(round(self.params[i][3],2)) + " \\textpm " +  str(self.errors[i][3]))
                    
                    try:
                        if self.params[i][6] == False: #Voigt not Double Lorentzian
                            Br_vals.append(str(round(self.params[i][4],2)) + " \\textpm " +  str(self.errors[i][4]))
                            SOS_vals.append(str(round(self.params[i][5],2)) + " \\textpm " +  str(self.errors[i][5]))
                            doublet = True

                    except:
                        pass                    
                        
                    try:    
                        if self.params[i][7] == False: #Double Lorentzian, not CK
                        
                            asym_vals.append(str(round(self.params[i][4],2)) + " \\textpm " +  str(self.errors[i][4]))
                            Br_vals.append(str(round(self.params[i][5],2)) + " \\textpm " +  str(self.errors[i][5]))
                            SOS_vals.append(str(round(self.params[i][6],2)) + " \\textpm " +  str(self.errors[i][6]))
                            doublet = True
                            asym = True
                    except:
                        pass
                    
                    try:    
                        if self.params[i][8] == True: #CK is fit as two singlets 
                          
                            asym_vals.append(str(round(self.params[i][4],2)) + " \\textpm " +  str(self.errors[i][4]))
                            Br_vals.append(str(round(self.params[i][5],2)) + " \\textpm " +  str(self.errors[i][5]))
                            SOS_vals.append(str(round(self.params[i][6],2)) + " \\textpm " +  str(self.errors[i][6]))
                            CK_vals.append(str(round(self.params[i][7],3)))
                            doublet = True
                            asym = True
                            CK = True
                    except:
                        pass
                    
                    try:
                        if self.params[i][4] >= 1: #If doublet and not DL this could be the Br instead --> Br should never be greater than 1 but asym should always be at least one or greater
                            if is_singlet == True: #CK
                                pass
                            elif is_singlet == False: #DL doublet
                                pass
                            else:
                          
                                asym_vals.append(str(round(self.params[i][4],2)) + " \\textpm " +  str(self.errors[i][4]))
                                asym = True
                    except:
                        pass 
                    area_vals.append(str(round(self.peak_areas[i],2)) + " \\textpm " +  str(round(self.upper_error_area[i],3)))
                    fwhm_vals.append(str(round(self.FWHM_values[i],2)))




            data["BE"] = BE_vals
            data["GFWHM"] = sigma_vals
            data["LFWHM"] = gamma_vals
            data["Amp"] = amp_vals
            
            if doublet == True:
                if asym == True:
                    data["asym"] = asym_vals
                data["Br"] = Br_vals
                data["SOS"] = SOS_vals
                if CK == True:
                    data["LFWHM_CK"] = CK_vals
            else: 
                if asym == True: #Double Lorentzian Singlet
                    data["asym"] = asym_vals


            data["Area"] = area_vals
            data["FWHM"] = fwhm_vals


            textabular = f"{'c'*len(headers)}"
            texheader = " & ".join(headers) + "\\\\"
            texdata = "\\hline\n"
            for label in data:
                texdata += f"{label} & {' & '.join(map(str,data[label]))} \\\\\n"



            table_begin = str("\\begin{tabular}{"+textabular+"}")
            table_header = str(texheader)
            hline = str("\\hline\n")
            table_data = str(texdata)
            table_end = str("\\end{tabular}")
            table_file.write(table_begin)
            table_file.write(hline)
            table_file.write(hline)
            table_file.write(table_header)
            table_file.write(table_data)
            table_file.write(hline)
            table_file.write(hline)
            table_file.write(table_end)
            table_file.close()
            print("Fitting Parameters Saved as LaTeX Table")

        def updateExportDataRows(args):
            #self.num = int(self.number_of_peaks.get())
            self.num = 1
            #self.Peak_number = abs(self.peak_number)

            i=int(self.peak_number.get()-1)
            #voigt_for_peak_type = 'Voigt'
            #values = self.peak_types in self.peakTypes_entries[i] if you want multiple background types


            #for row in range(1,((4*int(self.peak_number.get()))+1),4): #4 because 4 paramters. add changes here to output more values



            peakType = self.params[i][-1] #Last element in array is the curve fit type
            
           

            if(peakType.lower() == "lorentzian"):
                BE_text = str(round(self.params[i][0],2)) + "  " + "+/-" + "  " + str(self.errors[i][0])
                self.BE_text = StringVar(value=BE_text)

                lorentzian_text = str(round(self.params[i][1],3)) + "  " + "+/-" + "  " + str(self.errors[i][1])
                self.lorentzian_text = StringVar(value=lorentzian_text)

                amp_text = str(round(self.params[i][2],2)) + "  " + "+/-" + "  " + str(self.errors[i][2])
                self.amp_text = StringVar(value=amp_text)

                area_text = str(round(self.peak_areas[i],2)) + "  " + "+/-" + "  " + str(round(self.upper_error_area[i],3))
                self.area_text = StringVar(value=area_text)

                FWHM_text = str(round(self.FWHM_values[i],2))
                self.FWHM_text = StringVar(value=FWHM_text)

                self.peak_energy_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.BE_text, font=self.entryFont)
                self.peak_energy_entries[i].grid(column=1, row=2, sticky=(W, E))

                self.sigma_entries[i] = ttk.Label(self.analysis_tab, textvariable="0", font=self.entryFont)
                self.sigma_entries[i].grid(column=1, row=3, sticky=(W, E))

                self.lorentzian_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.lorentzian_text, font=self.entryFont)
                self.lorentzian_entries[i].grid(column=1, row=4, sticky=(W, E))

                self.amp_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.amp_text, font=self.entryFont)
                self.amp_entries[i].grid(column=1, row=5, sticky=(W, E))

                self.area_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.area_text, font=self.entryFont)
                self.area_entries[i].grid(column=1, row=6, sticky=(W, E))

                self.FWHM_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.FWHM_text, font=self.entryFont)
                self.FWHM_entries[i].grid(column=1, row=7, sticky=(W, E))

            elif(peakType.lower() == "gaussian"):
                BE_text = str(round(self.params[i][0],2)) + "  " + "+/-" + "  " + str(self.errors[i][0])
                self.BE_text = StringVar(value=BE_text)

                sigma_text = str(round(self.params[i][1],3)) + "  " + "+/-" + "  " + str(self.errors[i][1])
                self.sigma_text = StringVar(value=sigma_text)

                amp_text = str(round(self.params[i][2],2)) + "  " + "+/-" + "  " + str(self.errors[i][2])
                self.amp_text = StringVar(value=amp_text)

                area_text = str(round(self.peak_areas[i],2)) + "  " + "+/-" + "  " + str(round(self.upper_error_area[i],3))
                self.area_text = StringVar(value=area_text)

                FWHM_text = str(round(self.FWHM_values[i],2))
                self.FWHM_text = StringVar(value=FWHM_text)

                self.peak_energy_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.BE_text, font=self.entryFont)
                self.peak_energy_entries[i].grid(column=1, row=2, sticky=(W, E))

                self.sigma_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.sigma_text, font=self.entryFont)
                self.sigma_entries[i].grid(column=1, row=3, sticky=(W, E))

                self.lorentzian_entries[i] = ttk.Label(self.analysis_tab, textvariable="0", font=self.entryFont)
                self.lorentzian_entries[i].grid(column=1, row=4, sticky=(W, E))

                self.amp_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.amp_text, font=self.entryFont)
                self.amp_entries[i].grid(column=1, row=5, sticky=(W, E))

                self.area_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.area_text, font=self.entryFont)
                self.area_entries[i].grid(column=1, row=6, sticky=(W, E))

                self.FWHM_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.FWHM_text, font=self.entryFont)
                self.FWHM_entries[i].grid(column=1, row=7, sticky=(W, E))
            else: #Voigt, Double-Lorentzian, Doniach-Sunjic --> Asymmetry values not being shown right now


                BE_text = str(round(self.params[i][0],2)) + "  " + "+/-" + "  " + str(self.errors[i][0])
                self.BE_text = StringVar(value=BE_text)

                sigma_text = str(round(self.params[i][1],3)) + "  " + "+/-" + "  " + str(self.errors[i][1])
                self.sigma_text = StringVar(value=sigma_text)

                lorentzian_text = str(round(self.params[i][2],3)) + "  " + "+/-" + "  " + str(self.errors[i][2])
                self.lorentzian_text = StringVar(value=lorentzian_text)

                amp_text = str(round(self.params[i][3],2)) + "  " + "+/-" + "  " + str(self.errors[i][3])
                self.amp_text = StringVar(value=amp_text)

                area_text = str(round(self.peak_areas[i],2)) + "  " + "+/-" + "  " + str(round(self.upper_error_area[i],3))
                self.area_text = StringVar(value=area_text)

                FWHM_text = str(round(self.FWHM_values[i],2))
                self.FWHM_text = StringVar(value=FWHM_text)

                self.peak_energy_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.BE_text, font=self.entryFont)
                self.peak_energy_entries[i].grid(column=1, row=2, sticky=(W, E))

                self.sigma_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.sigma_text, font=self.entryFont)
                self.sigma_entries[i].grid(column=1, row=3, sticky=(W, E))

                self.lorentzian_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.lorentzian_text, font=self.entryFont)
                self.lorentzian_entries[i].grid(column=1, row=4, sticky=(W, E))

                self.amp_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.amp_text, font=self.entryFont)
                self.amp_entries[i].grid(column=1, row=5, sticky=(W, E))

                self.area_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.area_text, font=self.entryFont)
                self.area_entries[i].grid(column=1, row=6, sticky=(W, E))

                self.FWHM_entries[i] = ttk.Label(self.analysis_tab, textvariable=self.FWHM_text, font=self.entryFont)
                self.FWHM_entries[i].grid(column=1, row=7, sticky=(W, E))


            #i += 1

            #print(peak_labels)
        '''
        def get_areas(self,dir,params):
            #dir = str(dir.get())
            self.xps_analysis = xps_analysis2.xps_analysis(dir,params)
            self.xps_analysis.extract_data(plot_err=False)
            totalArea, peakAreas = self.xps_analysis.analyze()
            print(totalArea)
        '''
        """
        TODO:
        read in output file and get the relevant values to put in the boxes
        """
        self.analysis_obj = Analysis_plot(self.analysis_tab)

        # analysis_chisqr = StringVar(self.analysis_tab, 0.0)
        # analysis_background = StringVar(self.analysis_tab, 0.0)

        # For now put in placeholders
        # Select Folder Button
        analysis_button = ttk.Button(self.analysis_tab, text="Select Folder",command=select_analysis_folder)  # Add command to export data
        analysis_button.grid(column=0, row=0, sticky=(W, E), padx=self.padx, pady=self.pady, columnspan=2)

        self.peak_number = DoubleVar(self.root, 1) #Changed from 0 to 1 --> works for 1 peak not multiple. Combobox goes away. I think this is being used after folder is selected instead of waiting for user selection
        peak_number_options = [1,2,3,4,5,6,7,8,9,10]
        peak_number_entry = ttk.Combobox(self.analysis_tab, textvariable=self.peak_number, font=self.entryFont,values= peak_number_options)
        peak_number_entry.grid(column=1, row=1, sticky=(W, E)) #on same row as background checkbox
        peak_number_entry.bind('<<ComboboxSelected>>', updateExportDataRows)



        peak_labels = "Peak"
        BE_labels = "Energy"
        sigma_labels = "GFWHM"
        lorentzian_labels = "LFWHM"
        amp_labels = "Amp"
        area_labels = "Area"
        FWHM_labels = "FWHM"
        self.num_peak_labels = ttk.Label(self.analysis_tab, text=peak_labels, font=self.labelFont)
        self.num_peak_labels.grid_configure(column=0, row=1, sticky=W, padx=self.padx, pady=self.pady)
        self.BE_labels = ttk.Label(self.analysis_tab, text=BE_labels, font=self.labelFont)
        self.BE_labels.grid_configure(column=0, row=2, sticky=W, padx=self.padx, pady=self.pady)
        self.sigma_labels = ttk.Label(self.analysis_tab, text=sigma_labels, font=self.labelFont)
        self.sigma_labels.grid_configure(column=0, row=3, sticky=W, padx=self.padx, pady=self.pady)
        self.lorentzian_labels = ttk.Label(self.analysis_tab, text=lorentzian_labels, font=self.labelFont)
        self.lorentzian_labels.grid_configure(column=0, row=4, sticky=W, padx=self.padx, pady=self.pady)
        self.amp_labels = ttk.Label(self.analysis_tab, text=amp_labels, font=self.labelFont)
        self.amp_labels.grid_configure(column=0, row=5, sticky=W, padx=self.padx, pady=self.pady)
        self.area_labels = ttk.Label(self.analysis_tab, text=area_labels, font=self.labelFont)
        self.area_labels.grid_configure(column=0, row=6, sticky=W, padx=self.padx, pady=self.pady)
        self.FWHM_labels = ttk.Label(self.analysis_tab, text=FWHM_labels, font=self.labelFont)
        self.FWHM_labels.grid_configure(column=0, row=7, sticky=W, padx=self.padx, pady=self.pady)




        # Entries___________________________________________________________________________________________
        # entry_chisqr_best = ttk.Label(self.analysis_tab, textvariable=analysis_chisqr, font=self.entryFont, borderwidth=2,
        #                           relief="groove", background='#a9a9a9')
        # entry_chisqr_best.grid(column=1, row=4, sticky=(W, E), padx=self.padx)

        # entry_background = ttk.Label(self.analysis_tab, textvariable=analysis_background, font=self.entryFont,borderwidth=2,
        #                               relief="groove", background='#a9a9a9')
        # entry_background.grid(column=1, row=5, sticky=(W, E), padx=self.padx)

        # Number of Peak to Find Row
        self.numPeaks = int(self.number_of_peaks.get())

        # Chisqr Entry
        # entry_chisqr_best = ttk.Label(self.analysis_tab, textvariable=analysis_chisqr, font=self.entryFont, borderwidth=2,
        #                           relief="groove", background='#a9a9a9')
        # entry_chisqr_best.grid(column=1, row=4, sticky=(W, E), padx=self.padx)

        # Plot Best Fit Button
        button_plot = ttk.Button(self.analysis_tab,text="Plot Best Fit",command=calculate_and_plot)  # Add command to plot data using postprocessing
        button_plot.grid(column=0, row=((4*self.numPeaks)+5), sticky=(W, E), padx=self.padx, pady=self.pady, columnspan=2)



        # Export Values Button

        button_export = ttk.Button(self.analysis_tab, text="Export Values to LaTeX Table", command = export_peak_parameters)  # Add command to export data
        button_export.grid(column=0, row=((4*self.numPeaks)+6), sticky=(W, E), padx=self.padx, pady=self.pady, columnspan=2)

        button_export_fit = ttk.Button(self.analysis_tab, text="Export Peak Fit", command = export_indi_peak_fit)  # Add command to export data
        button_export_fit.grid(column=0, row=((4*self.numPeaks)+7), sticky=(W, E), padx=self.padx, pady=self.pady, columnspan=2)

        self.analysis_tab.columnconfigure(3, weight=1)
        self.analysis_tab.rowconfigure(0, weight=1)
        self.analysis_tab.rowconfigure(1, weight=1)
        self.analysis_tab.rowconfigure(2, weight=1)
        self.analysis_tab.rowconfigure(3, weight=1)
        self.analysis_tab.rowconfigure(4, weight=1)
        self.analysis_tab.rowconfigure(5, weight=1)
        self.analysis_tab.rowconfigure(6, weight=1)
        self.analysis_tab.rowconfigure(7, weight=1)
        self.analysis_tab.rowconfigure(8, weight=1)
        self.analysis_tab.rowconfigure(9, weight=1)

    def stop_term(self):
        #print("In stop term")
        #print("PID TO KILL ", self.pid_list[len(self.pid_list)-1])
        #command = 'kill -9 ', self.pid_list[len(self.pid_list)-1]
        #os.killpg(self.pid_list[len(self.pid_list)-1], signal.SIGTERM)
        #subprocess.Popen(command)
        # print("/n/n/n/n/n/n/n/n/n/n"
        #       "**************"
        #       "*******************"
        #       "(++++&****************************"
        #       ""
        #       "/n/n/n/n/n/n")

        self.stop_not_pressed = False
        if not self.command_list.empty():
            self.stop_all()
        elif hasattr(self, 'proc'):
            self.proc.kill()
            self.proc.terminate()
            print("Stopped xps_neo")
            self.proc.kill()

    def on_closing(self):
        """
        on closing function
        """
        if messagebox.askokcancel("Quit", "Do you want to quit?"):
            self.stop_term()
            if hasattr(self, 'terminal'):
                self.root.quit()
                self.terminal.destroy()
            else:
                self.root.quit()

    def Run(self):
        """
        Run the code
        """
        self.root.protocol('WM_DELETE_WINDOW', self.on_closing)
        self.root.mainloop()


def main():
    """Console entry point (pyproject: xps_neo_gui)."""
    root = App()
    root.Run()


if __name__ == "__main__":
    main()
