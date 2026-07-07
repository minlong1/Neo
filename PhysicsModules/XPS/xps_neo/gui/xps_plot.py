import tkinter as tk

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
#import NanoIndent_Analysis
from PhysicsModules.XPS.xps_neo.gui import xps_analysis2

class Analysis_plot:
    def __init__(self, frame, data_KE = False, data_XES = False):
        
        self.frame = frame
        self.fig = Figure(figsize=(6.4,4.8), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame)
        # Create initial figure canvas
        self.canvas.get_tk_widget().grid(column=2, row=1, rowspan=13, columnspan=5, sticky="nsew",
                                         padx=5, pady=5)
        self.ax = self.fig.add_subplot(111)
        # create toolbar
        self.toolbarFrame = tk.Frame(master=self.frame)
        self.toolbarFrame.grid(column=2, row=0, rowspan=1, columnspan=5, sticky="nsew")
        toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
        self.params = {}

    def unique_parameters(self,dir,params,scale_var,data_KE,data_XES,peakType, title):
        dir = str(dir.get())
        self.xps_analysis = xps_analysis2.xps_analysis(dir,params, peakType)
        self.xps_analysis.extract_data(scale_var,plot_err=False) #Maybe add scale_var to this function to make it seen in xps_individual?
        self.xps_analysis.get_unique(self)

    def initial_parameters(self,dir,params,scale_var,data_KE,data_XES,peakType, title):
        dir = str(dir.get())
        self.xps_analysis = xps_analysis2.xps_analysis(dir,params, peakType)
        self.xps_analysis.extract_data(scale_var,plot_err=False) #Maybe add scale_var to this function to make it seen in xps_individual?
        self.xps_analysis.score()
        self.totalArea, self.areas = self.xps_analysis.analyze()
        self.xps_analysis.analyze()
        self.fig.clf()
        self.xps_analysis.plot_data(scale_var, data_KE, data_XES, title=title,fig_gui = self.fig)
        self.canvas.draw()
        return self.xps_analysis.get_params()
        '''
        self.nano_analysis.score()
        self.nano_analysis.calculate_parameters(verbose=False)
        self.fig.clf()
        self.nano_analysis.plot_data(title=title,fig_gui=self.fig)
        self.canvas.draw()
        return self.nano_analysis.get_params()
        '''
class Data_plot:
    def __init__(self, frame, data_KE = False, data_XES = False):
        self.frame = frame
        self.fig = Figure(figsize=(6.4, 4.8), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame)
        # Create initial figure canvas
        self.canvas.get_tk_widget().grid(column=0, row=1, rowspan=4, columnspan=8, sticky="nsew",
                                         padx=5, pady=5)
        self.ax = self.fig.add_subplot(111)
        # create toolbar
        self.toolbarFrame = tk.Frame(master=self.frame)
        self.toolbarFrame.grid(column=3, row=0, rowspan=1, columnspan=5, sticky="w")
        toolbar = NavigationToolbar2Tk(self.canvas, self.toolbarFrame)
        self.params = {}

    def initial_parameters(self,data_obj, data_KE, data_XES):
        """
        fileraw
        """
        self.data_obj = data_obj
        self.title = 0
        
        if data_KE == False:
            self.xlabel = 'Binding Energy (eV)'
     
        else:
            self.xlabel = 'Kinetic Energy (eV)'

        if data_XES == True:
            self.xlabel = 'Photon Energy (eV)'
       
        self.ylabel= 'Intensity (a.u.)'
    #def plot written by evan 
    def plot(self,x_data_array,y_data_array, param_label,title, scale_var,data_KE, data_XES):
        x = self.data_obj.get_x(data_KE, data_XES)
        y = self.data_obj.get_y(scale_var)


        '''
        min_y = np.min(y)
       
        for i in range(len(y)):
            y[i] = y[i] - min_y
        '''
        
        self.ax.clear()
        '''
        for i in range(len(x_data_array)):
            
            self.ax.plot(x_data_array[i],y_data_array[i], 'b.', label = param_label)
        '''
        #global data_KE
        self.ax.plot(x,y,color='gray', linestyle='None', marker='o', markeredgewidth='0.3', markerfacecolor='None', markeredgecolor='gray', markersize='5', label='Data')
        
        if data_KE == False:
            self.ax.invert_xaxis()   
        else:
            pass
        if data_XES == True:
            self.ax.invert_xaxis()  
        

        #self.ax.invert_xaxis()
        self.ax.legend()
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        self.ax.set_title(title)
        self.fig.tight_layout()

        self.canvas.draw()

    