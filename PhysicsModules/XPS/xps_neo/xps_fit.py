"""
Author: Evan Restuccia (evan@restuccias.com)
"""
import numpy as np
import sys
import scipy as scipy
from scipy import signal
import random
import decimal
import matplotlib.pyplot as plt
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning #Use this to get rid of numba warnings that are outputted in the terminal during run
import warnings
import numba as nb
import pathlib
import re



warnings.filterwarnings('ignore')
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)
#import chart_studio.plotly as py
#import peakutils
#from pybaselines import Baseline, utils #Will have to have people install this in order to make it run

class peak():
    """One fittable XPS peak (Voigt/Gaussian/Lorentzian/Double Lorentzian/
    Doniach-Sunjic; singlet, spin-orbit doublet, or Coster-Kronig), holding
    its parameter ranges/guesses and the shape function selected by
    `peakType` (`self.peakFunc`). See `__init__`'s docstring for the
    free-vs-semi-free/correlated parameter model and how to add a new
    peak type."""

    #Probably default these to -1 later as a test condition, but not yet sure
    def __init__(self,paramRange,peakType,is_singlet = 'is_singlet', is_coster_kronig = 'is_coster_kronig', BE_correlated = "BE_correlated", BE_correlated_mult = "BE_correlated_mult", sigma_correlated = "sigma_correlated", sigma_correlated_mult = "sigma_correlated_mult", gamma_correlated = "gamma_correlated", gamma_correlated_mult = "gamma_correlated_mult", amp_correlated = "amp_correlated", amp_correlated_mult = "amp_correlated_mult", asym_correlated = "asym_correlated", asym_correlated_mult = "asym_correlated_mult", sos_correlated = "sos_correlated", sos_correlated_mult = "sos_correlated_mult", br_correlated = "br_correlated", br_correlated_mult = "br_correlated_mult"):

        """
        takes in dictionary paramRange, with names of parameters and their corresponding allowed values
        takes the dictionary and sets the ranges for each
        Two types of params: Free and semi-free
        Free params get fully random values within their range
        Semi-Free take in a starting guess and modify it within an allowed range

        Then initalizes the peakType and hooks the proper function up to the correspond peakFunc

        To add a Peak type:
        Add its peakType option to the picker in the GUI
        Give its parameters ranges and pass them into the write_ini in the gui
        add the range to ini_parser and also the paramrange dict in xps_neo
        add that bkgnType in init here, and assign the function to self.peakFunc
        add the get function for it
        add the set function in the gui xps_fit
        """
        #fetch ranges for the values from the dict
        self.paramRange= paramRange
        self.is_singlet = is_singlet
        self.is_coster_kronig = is_coster_kronig
        self.BE_correlated = BE_correlated
        self.BE_correlated_mult = BE_correlated_mult
        self.sigma_correlated = sigma_correlated
        self.sigma_correlated_mult = sigma_correlated_mult
        self.gamma_correlated = gamma_correlated
        self.gamma_correlated_mult = gamma_correlated_mult
        self.amp_correlated = amp_correlated
        self.amp_correlated_mult = amp_correlated_mult
        self.asym_correlated = asym_correlated
        self.asym_correlated_mult = asym_correlated_mult
        self.sos_correlated = sos_correlated
        self.sos_correlated_mult = sos_correlated_mult
        self.br_correlated = br_correlated
        self.br_correlated_mult = br_correlated_mult

        if paramRange == '':
            # Analysis container mode (Phase 3d): parameters are assigned
            # explicitly via set()/set_* afterwards; no ranges are built and
            # no RNG is consumed. Numeric defaults are placeholders only.
            self.bindingEnergy = 0.0
            self.gaussian = 0.0
            self.lorentz = 0.0
            self.amp = 0.0
            self.asymmetry = 1.0
            self.asymmetryDoniach = 0.0
            self.branching_ratio = 0.5
            self.spinOrbitSplit = 0.0
            self.lorentz_CK = 0.0
            self.asym_CK = 0.0
            self.peak_y = []
            self.peakType = peakType
            if(self.peakType.lower() == "voigt"):
                self.func = self.voigtFunc
            elif(self.peakType.lower() == "gaussian"):
                self.func = self.gaussFunc
            elif(self.peakType.lower() == "lorentzian"):
                self.func = self.lorentzFunc
            elif(self.peakType.lower() == "double lorentzian"):
                self.func = self.doubleLorentzFunc
            elif(self.peakType.lower() == "doniach-sunjic"):
                self.func = self.doniachSunjicFunc
            else:
                print("Error assigning peak type")
                print("Peaktype found is: " + str(self.peakType))
                sys.exit(1)
            return

        try:
            self.gaussRange = np.arange(paramRange['Gaussian'][0],paramRange['Gaussian'][1],paramRange['Gaussian'][2])

            self.lorentzRange = np.arange(paramRange['Lorentzian'][0],paramRange['Lorentzian'][1],paramRange['Lorentzian'][2])

            self.lorentzCKRange = np.arange(paramRange['Lorentzian Coster-Kronig'][0],paramRange['Lorentzian Coster-Kronig'][1],paramRange['Lorentzian Coster-Kronig'][2])
            
            self.bindingEnergyRange = np.arange(paramRange['Binding Energy'][0],paramRange['Binding Energy'][1],paramRange['Binding Energy'][2])
            #print("FIT", paramRange['Amplitude'])
            self.ampRange = np.arange(paramRange['Amplitude'][0],paramRange['Amplitude'][1],paramRange['Amplitude'][2])
            self.asymmetryRange = np.arange(paramRange['Asymmetry'][0],paramRange['Asymmetry'][1],paramRange['Asymmetry'][2])

            self.asymCKRange = np.arange(paramRange['Asymmetry Range Coster-Kronig'][0],paramRange['Asymmetry Range Coster-Kronig'][1],paramRange['Asymmetry Range Coster-Kronig'][2])

            self.asymmetryDoniachRange = np.arange(paramRange['Asymmetry Doniach-Sunjic'][0],paramRange['Asymmetry Doniach-Sunjic'][1],paramRange['Asymmetry Doniach-Sunjic'][2])
       
            try:
                self.spinOrbitSplitRange = np.arange(paramRange['spinOrbitSplitting'][0],paramRange['spinOrbitSplitting'][1],paramRange['spinOrbitSplitting'][2])
            except:
                self.spinOrbitSplitRange = [0,0,0]
            try:
                self.branchingRatioRange = np.arange(paramRange['br_range'][0],paramRange['br_range'][1],paramRange['br_range'][2])
            except:
                self.branchingRatioRange = [0,0,0]

            #fully free within their range
            self.gaussian = np.random.choice(self.gaussRange)
            self.lorentz = np.random.choice(self.lorentzRange)
            self.amp = np.random.choice(self.ampRange)
            self.asymmetry = np.random.choice(self.asymmetryRange)
            self.asymmetryDoniach = np.random.choice(self.asymmetryDoniachRange)
          
            
            #IF CORRELATED MAKE THE NP.RANDOM.CHOICE BE THE SAME AS THE CORRELATED PEAK??? ---HOW

            #self.lorentz_CK = np.random.choice(self.lorentzCKRange) #Coster-Kronig
            
        except:
            print("Error in init of xps_fit")
            sys.exit(1)


        #the range is a modifier on the input value
        self.bindingEnergy= np.random.choice(self.bindingEnergyRange) #Not using rn --> changed to set parameter
      
     
        #self.s_o_split = np.random.choice(self.s_o_splittingRange) #Not using rn

        self.peakType = peakType
        self.peak_y = []
        

        #self.is_singlet = is_singlet
        #self.branching_ratio = paramRange['branching_ratio']
        self.branching_ratio = np.random.choice(self.branchingRatioRange)
        #self.lorentz_CK = self.lorentz/self.branching_ratio #Constraining the lorentz_CK to be proportional to the lorentz value based on proper branching ratio --> Should equalize areas



        #Coster-Kronig Values:

        self.lorentz_CK =   np.random.choice(self.lorentzCKRange) #Change range so that self.lorentz_CK is always bigger than self.lorentz
        
        if self.lorentz_CK < self.lorentz:
            self.lorentz_CK = self.lorentz + 0.01

        self.asym_CK =   np.random.choice(self.asymCKRange) #Change range so that self.lorentz_CK is always bigger than self.lorentz
        
        if self.asym_CK < self.asymmetry:
            self.asym_CK = self.asymmetry + 0.01
        
         


        #self.spinOrbitSplit = paramRange['spinOrbitSplit']
        #self.spinOrbitSplitRange = np.arange(paramRange['Spin-Orbit Split'][0],paramRange['Spin-Orbit Split'][1],paramRange['Spin-Orbit Split'][2])
        self.spinOrbitSplit= np.random.choice(self.spinOrbitSplitRange)
        self.peakType = peakType
        if(self.peakType.lower() == "voigt"):
            self.func = self.voigtFunc
        elif(self.peakType.lower() == "gaussian"):
            self.func = self.gaussFunc
        elif(self.peakType.lower() == "lorentzian"):
            self.func = self.lorentzFunc
        elif(self.peakType.lower() == "double lorentzian"):
            self.func = self.doubleLorentzFunc
        elif(self.peakType.lower() == "doniach-sunjic"):
            self.func = self.doniachSunjicFunc

        else:
            print("Error assigning peak type")
            print("Peaktype found is: " + str(self.peakType))
            sys.exit(1)





    # ---- Phase 4.2: slim pickling for --workers -------------------------
    # Workers only evaluate curves; the np.arange mutation ranges (the bulk
    # of the object) and the bound function pointer are dropped from the
    # pickle and the pointer is re-hooked on unpickle. Mutation always runs
    # in the main process, where the original arrays stay untouched
    # (_absorb_state never overwrites attributes absent from the worker copy).
    _PICKLE_SKIP = ("func", "gaussRange", "lorentzRange", "lorentzCKRange",
                    "bindingEnergyRange", "ampRange", "asymmetryRange",
                    "asymCKRange", "asymmetryDoniachRange",
                    "spinOrbitSplitRange", "branchingRatioRange")

    def __getstate__(self):
        return {k: v for k, v in self.__dict__.items()
                if k not in self._PICKLE_SKIP}

    def __setstate__(self, state):
        self.__dict__.update(state)
        if(self.peakType.lower() == "voigt"):
            self.func = self.voigtFunc
        elif(self.peakType.lower() == "gaussian"):
            self.func = self.gaussFunc
        elif(self.peakType.lower() == "lorentzian"):
            self.func = self.lorentzFunc
        elif(self.peakType.lower() == "double lorentzian"):
            self.func = self.doubleLorentzFunc
        elif(self.peakType.lower() == "doniach-sunjic"):
            self.func = self.doniachSunjicFunc

    def peakFunc(self,x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split):
        #print(self.func(x), "HELLO")
        return self.func(x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split)

    #----------------------------Getters-------------------------------------------------------



    def get(self,for_output_data_file = False):
        """
        Gets params in format [Singlet Params, bool for is_singlet, doublet params, peakType]
        """
        #Add in FWHM values for each curve type --> Voigt has 3, gauss and loretnz only have one
        params = []
        if self.peakType.lower() == 'voigt':
            params = [self.bindingEnergy,self.gaussian,self.lorentz,self.amp] #mutate relies on the order here, so to change this you need to change
            #params = [self.bindingEnergy,self.gaussian,self.lorentz,self.amp,self.peakType] #mutate relies on the order here, so to change this you need to change mutate
        if self.peakType.lower() == 'gaussian':
            params = [self.bindingEnergy,self.gaussian,self.amp]#everything except lorentzian which in this case is the width i believe --> Still want FWHM
            #params = [self.bindingEnergy,self.gaussian,self.amp,self.peakType]#everything except lorentzian which in this case is the width i believe --> Still want FWHM
        if self.peakType.lower() == 'lorentzian':
            params = [self.bindingEnergy,self.lorentz,self.amp]
            #params = [self.bindingEnergy,self.lorentz,self.amp,self.peakType]
        if self.peakType.lower() == 'double lorentzian':
            params = [self.bindingEnergy,self.gaussian,self.lorentz,self.amp,self.asymmetry]
            #params = [self.bindingEnergy,self.gaussian,self.lorentz,self.amp,self.asymmetry,self.peakType]
        if self.peakType.lower() == 'doniach-sunjic':
            params = [self.bindingEnergy,self.gaussian,self.lorentz,self.amp,self.asymmetryDoniach] #no gaussian needed
            #params = [self.bindingEnergy,self.lorentz,self.amp,self.asymmetry,self.peakType] #no gaussian needed

        #grab the bool is_singlet if writing for output
        if(for_output_data_file):
            params.append(self.is_singlet)
          

        #grab the doublet params if its a doublet
        if not(self.is_singlet):
            #if(for_output_data_file):
                #params.append(self.branching_ratio)
            params.append(self.branching_ratio)
            params.append(self.spinOrbitSplit)

            if self.is_coster_kronig: #If it is a doublet AND is Coster-Kronig --> I think it is not seeing CK correctly here... Throws a True value for if CK is selected... FIXlorentz_CK
                params.append(self.lorentz_CK)
                if self.peakType.lower() == 'double lorentzian':
                    params.append(self.asym_CK)
            
                #params.append(self.is_coster_kronig)
            params.append(self.is_coster_kronig)


            

        #and always end on the peakType
        params.append(self.peakType)


        #Add this peaks shirley background if using peak shirley
        if len(params) == 0:
            print("Cant do 'def get' in peaks class in XPS_FIT, most likely a new peak was added and needs to be added to the get options")
            sys.exit(1)
        else:
            return params



    """
    I don't think this set of functions is used ever
    """
    def getGaussian(self):
        return self.gaussian
    def getLorenztian(self):
        return self.lorentz
    def getAmplitude(self):
        return self.amp
    def getBindingEnergy(self):
        return self.bindingEnergy
    def getAsymmetry(self):
        return self.asymmetry
    def getAsymmetryDoniach(self):
        return self.asymmetryDoniach
    def getSpinOrbitSplit(self):
        return self.spinOrbitSplit
    def getBranchingRatio(self):
        return self.branching_ratio




    def getFWHM(self,x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split): #Not wure where to call this in order to output the FWHM values for each peak curve type
        self.bindingEnergy = BE
        self.lorentz = width
        self.lorentz_CK = width_CK
        self.is_coster_kronig = coster_kronig
        self.gaussian = sigma
        self.amp = A
        self.asymmetry = asym
        self.asymmetryDoniach = asymD
        self.is_singlet = singlet
        self.branching_ratio = branch
        self.spinOrbitSplit = split
       
        #MAD = np.mean(np.absolute(self.peak_y-np.mean(self.y))) #Mean Absolute Deviaition
        #self.lorentz = gamma = FWHM of a pure Lorentzian lineshape
        if(self.peakType.lower() == "voigt"):
            self.fwhm_g = self.gaussian #2*self.gaussian*np.sqrt(2*np.log(2))
            self.fwhm_l = self.lorentz
            self.fwhm_v = 0.5346*self.fwhm_l + np.sqrt(0.2166*pow(self.fwhm_l,2) + pow(self.fwhm_g, 2))
            #self.fwhm_v = self.fwhm_l/2 + np.sqrt(pow(self.fwhm_l, 2)/4 + pow(self.fwhm_g, 2))
            return self.fwhm_v #, self.fwhm_l, self.fwhm_g
        elif(self.peakType.lower() == "gaussian"):
            self.fwhm_g = self.gaussian #2*self.gaussian*np.sqrt(2*np.log(2))
            #print(self.fwhm_g)
            return self.fwhm_g
        elif(self.peakType.lower() == "lorentzian"):
            self.fwhm_l = self.lorentz
            return self.fwhm_l
        elif(self.peakType.lower() == "double lorentzian"):

            if self.is_coster_kronig == True: 
                self.fwhm_g = self.gaussian #2*self.gaussian*np.sqrt(2*np.log(2))
                self.fwhm_l_CK = self.lorentz_CK
                fwhm_asym_CK = self.lorentz_CK*self.asymmetry
                self.fwhm_1_CK = 0.5346*self.fwhm_l_CK + np.sqrt(0.2166*pow(self.fwhm_l_CK,2) + pow(self.fwhm_g, 2))
                self.fwhm_2_CK = 0.5346*fwhm_asym_CK + np.sqrt(0.2166*pow(fwhm_asym_CK,2) + pow(self.fwhm_g, 2))
                self.fwhm_dl_CK = self.fwhm_1_CK/2 + self.fwhm_2_CK/2 #Adding hwhm's together to get fwhm for total peak (including asymmetry) 


            self.fwhm_g = self.gaussian #2*self.gaussian*np.sqrt(2*np.log(2))
            self.fwhm_l = self.lorentz
            fwhm_asym = self.lorentz*self.asymmetry
            self.fwhm_1 = 0.5346*self.fwhm_l + np.sqrt(0.2166*pow(self.fwhm_l,2) + pow(self.fwhm_g, 2))
            self.fwhm_2 = 0.5346*fwhm_asym + np.sqrt(0.2166*pow(fwhm_asym,2) + pow(self.fwhm_g, 2))
            self.fwhm_dl = self.fwhm_1/2 + self.fwhm_2/2 #Adding hwhm's together to get fwhm for total peak (including asymmetry) 
            return self.fwhm_dl #, self.fwhm_l, self.fwhm_g #Should double lorentzian still provide FWHM_L?
        elif(self.peakType.lower() == "doniach-sunjic"):
            self.fwhm_g = self.gaussian #2*self.gaussian*np.sqrt(2*np.log(2))
            self.fwhm_l = self.lorentz #Should this be two values? One for the asymmetrical side of the peak
            self.fwhm_v = 0.5346*self.fwhm_l + np.sqrt(0.2166*pow(self.fwhm_l,2) + pow(self.fwhm_g, 2))
            #self.fhwm_v = self.fwhm_l/2 + np.sqrt(pow(self.fwhm_l, 2)/4 + pow(self.fwhm_g, 2))
            return self.fwhm_v #, self.fwhm_l, self.fwhm_g #Should double lorentzian still provide FWHM_G?
        else:
            print("Error in FWHM caluclation")
            print("Peaktype found is: " + str(self.peakType))
            pass

    def getY(self,x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split):
        """
        Gets the y values for a peak by running the peakFunc
        """
        #self.FWHM_values = self.getFWHM()
        #print("FWHM: ",FWHM_values)

        self.peakFunc(x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split)
      
        return self.peak_y

     #------------------------Setters--------------------------------------

    def set(self,params,doublet_params, CK_params, is_singlet=False, is_coster_kronig=False):
        """
        This is the main differentiator between the fit in gui and xpsfolders
        The set function takes in a list of params, and should expect them to come in the same order that they were pushed out
        in the get function.

        """
        if(self.peakType.lower() == "voigt"):
            self.set_voigt(params)
        elif(self.peakType.lower() == "gaussian"):
            self.set_gauss(params)
        elif(self.peakType.lower() == "lorentzian"):
            self.set_lorentz(params)
        elif(self.peakType.lower() == "double lorentzian"):
            self.set_doubleLorentz(params)
        elif(self.peakType.lower() == "doniach-sunjic"):
            self.set_doniachSunjic(params)
        else:
            print("Error, cant find type to set, type found was: " + str(self.peakType))

        # Doublet params apply to doublets (gui-fork semantics adopted in
        # Phase 3d; the pre-unification `if(is_singlet)` read doublet params
        # for singlets and crashed - set() is only called by the analysis).
        if not is_singlet:
            self.is_singlet = is_singlet
            self.branching_ratio = doublet_params[0]
            self.spinOrbitSplit = doublet_params[1]
            if is_coster_kronig == True:
                self.lorentz_CK = CK_params[0]
                if self.peakType.lower() == "double lorentzian":
                    self.asym_CK = CK_params[1]

    #probably can delete these function
    def setGaussian(self,newVal):
        self.gaussian = newVal
    def setLorentzian(self,newVal):
        self.lorentz = newVal
    def setAmplitude(self,newVal):
        self.amp = newVal
    def setBindingEnergy(self,newVal):
        self.bindingEnergy = newVal
    def setAsymmetry(self,newVal):
        self.asymmetry = newVal
    def setAsymmetryDoniach(self,newVal):
        self.asymmetryDoniach = newVal
    def setSpinOrbitSplit(self,newVal):
        self.spinOrbitSplit = newVal




    def set_voigt(self,paramList):
        self.bindingEnergy = paramList[0]
        self.gaussian = paramList[1]
        self.lorentz = paramList[2]
        self.amp = paramList[3]
        #self.FWHM_values = paramList[4]

    def set_gauss(self,paramList):
        self.bindingEnergy = paramList[0]
        self.gaussian = paramList[1]
        self.amp = paramList[2]
        #self.FWHM_values = paramList[3]
        #does lorentz effect gauss? --> it is just the width which is not used in the equation at all

    def set_lorentz(self,paramList):
        self.bindingEnergy = paramList[0]
        self.lorentz = paramList[1]
        self.amp = paramList[2]
        #self.FWHM_values = paramList[3]

    def set_doubleLorentz(self,paramList):
        self.bindingEnergy = paramList[0]
        self.gaussian = paramList[1]
        self.lorentz = paramList[2]
        self.amp = paramList[3]
        self.asymmetry = paramList[4]
        #self.FWHM_values = paramList[5]

    def set_doniachSunjic(self,paramList):
        self.bindingEnergy = paramList[0]
        self.gaussian = paramList[1]
        self.lorentz = paramList[2]
        self.amp = paramList[3]
        self.asymmetryDoniach = paramList[4]
        #self.FWHM_values = paramList[5]


    def checkOutbound(self, peak_num):
        """Check if out of bounds
        """
        #SOS, BR, and Coster-Kronig ranges??
        if self.peakType.lower() == 'voigt':

            BE_guess = self.paramRange['BE'][peak_num]
            amp_guess = self.paramRange['Amp'][peak_num]
            gamma_guess = self.paramRange['Gamma'][peak_num]
            sigma_guess = self.paramRange['Sigma'][peak_num]


            self.bindingEnergy = np.clip(self.bindingEnergy,BE_guess + self.paramRange['Binding Energy'][0],BE_guess + self.paramRange['Binding Energy'][1])
            self.gaussian = np.clip(self.gaussian,sigma_guess + self.paramRange['Gaussian'][0],sigma_guess + self.paramRange['Gaussian'][1])
            self.lorentz = np.clip(self.lorentz,gamma_guess + self.paramRange['Lorentzian'][0],gamma_guess + self.paramRange['Lorentzian'][1])
            self.amp = np.clip(self.amp,amp_guess + self.paramRange['Amplitude'][0],amp_guess + self.paramRange['Amplitude'][1])


        elif self.peakType.lower() == 'gaussian':
            BE_guess = self.paramRange['BE'][peak_num]
            amp_guess = self.paramRange['Amp'][peak_num]
            sigma_guess = self.paramRange['Sigma'][peak_num]

            self.bindingEnergy = np.clip(self.bindingEnergy,BE_guess + self.paramRange['Binding Energy'][0],BE_guess + self.paramRange['Binding Energy'][1])
            self.gaussian = np.clip(self.gaussian,sigma_guess + self.paramRange['Gaussian'][0],sigma_guess + self.paramRange['Gaussian'][1])
            self.amp = np.clip(self.amp,amp_guess + self.paramRange['Amplitude'][0],amp_guess + self.paramRange['Amplitude'][1])


        elif self.peakType.lower() == 'lorentzian':
            BE_guess = self.paramRange['BE'][peak_num]
            amp_guess = self.paramRange['Amp'][peak_num]
            gamma_guess = self.paramRange['Gamma'][peak_num]

            self.bindingEnergy = np.clip(self.bindingEnergy,BE_guess + self.paramRange['Binding Energy'][0],BE_guess + self.paramRange['Binding Energy'][1])
            self.lorentz = np.clip(self.lorentz,gamma_guess + self.paramRange['Lorentzian'][0],gamma_guess + self.paramRange['Lorentzian'][1])
            self.amp = np.clip(self.amp,amp_guess + self.paramRange['Amplitude'][0],amp_guess + self.paramRange['Amplitude'][1])

        elif self.peakType.lower() == 'double lorentzian':
            BE_guess = self.paramRange['BE'][peak_num]
            amp_guess = self.paramRange['Amp'][peak_num]
            gamma_guess = self.paramRange['Gamma'][peak_num]
            sigma_guess = self.paramRange['Sigma'][peak_num]
            asym_guess = self.paramRange['Asym'][peak_num]


            self.bindingEnergy = np.clip(self.bindingEnergy,BE_guess + self.paramRange['Binding Energy'][0],BE_guess + self.paramRange['Binding Energy'][1])
            self.gaussian = np.clip(self.gaussian,sigma_guess + self.paramRange['Gaussian'][0],sigma_guess + self.paramRange['Gaussian'][1])
            self.lorentz = np.clip(self.lorentz,gamma_guess + self.paramRange['Lorentzian'][0],gamma_guess + self.paramRange['Lorentzian'][1])
            self.amp = np.clip(self.amp,amp_guess + self.paramRange['Amplitude'][0],amp_guess + self.paramRange['Amplitude'][1])
            self.asymmetry = np.clip(self.asymmetry,asym_guess + self.paramRange['Asymmetry'][0],asym_guess + self.paramRange['Asymmetry'][1])
            
           

        elif self.peakType.lower() == 'doniach-sunjic':
            BE_guess = self.paramRange['BE'][peak_num]
            amp_guess = self.paramRange['Amp'][peak_num]
            gamma_guess = self.paramRange['Gamma'][peak_num]
            sigma_guess = self.paramRange['Sigma'][peak_num]

            self.bindingEnergy = np.clip(self.bindingEnergy,BE_guess + self.paramRange['Binding Energy'][0],BE_guess + self.paramRange['Binding Energy'][1])
            self.gaussian = np.clip(self.gaussian,sigma_guess + self.paramRange['Gaussian'][0],sigma_guess + self.paramRange['Gaussian'][1])
            self.lorentz = np.clip(self.lorentz,gamma_guess + self.paramRange['Lorentzian'][0],gamma_guess + self.paramRange['Lorentzian'][1])
            self.amp = np.clip(self.amp,amp_guess + self.paramRange['Amplitude'][0],amp_guess + self.paramRange['Amplitude'][1])
            self.asymmetry = np.clip(self.asymmetry,self.paramRange['Asymmetry Doniach-Sunjic'][0],self.paramRange['Asymmetry Doniach-Sunjic'][1])


    #---------------------------------------Mutation Functions-------------------------------------------------


    def mutate(self,chance):
        self.mutateGauss(chance)
        self.mutateAmplitude(chance)
        self.mutateBE(chance)
        self.mutateLorentz(chance)
        #self.mutateSplitting(chance)
        self.mutateAsymmetry(chance)
        self.mutateAsymmetryDoniach(chance)
        '''
        if is_singlet==False:
            self.mutateSpinOrbitSplit
        '''
    def mutateGauss(self,chance):
        if random.random()*100 < chance:
            self.gaussian = np.random.choice(self.gaussRange)
    def mutateLorentz(self,chance):
        if random.random()*100 < chance:
            self.lorentz = np.random.choice(self.lorentzRange)
    def mutateAmplitude(self,chance):
        if random.random()*100 < chance:
            self.amp = np.random.choice(self.ampRange)
    def mutateBE(self,chance):
        if random.random()*100 < chance:
            self.bindingEnergy = np.random.choice(self.bindingEnergyRange)
    def mutateSplitting(self,chance):
        if random.random()*100 < chance:
            self.s_o_split = np.random.choice(self.s_o_splittingRange)
    def mutateAsymmetry(self,chance):
        if random.random()*100 < chance:
            self.asymmetry = np.random.choice(self.asymmetryRange)
    def mutateAsymmetryDoniach(self,chance):
        if random.random()*100 < chance:
            self.asymmetryDoniach = np.random.choice(self.asymmetryDoniachRange)
    
    def mutateSpinOrbitSplit(self,chance):
        if random.random()*100 < chance:
            self.spinOrbitSplit = np.random.choice(self.spinOrbitSplitRange)

    def mutateBranchingRatio(self,chance):
        if random.random()*100 < chance:
            self.branching_ratio = np.random.choice(self.branchingRatioRange)
    
    #Peak curve fit equations start at line 6681 in Aanalyzer PUnit1


    #----------------------Peak Curve Form Definitions------------------------------------------#



    #A bit scrappy at the moment, may need cleaning later
   
    def voigtFunc(self,x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split):


        #Check to see if the data is in KE energy
        self.bindingEnergy = BE
        self.lorentz = width
        self.lorentz_CK = width_CK
        self.is_coster_kronig = coster_kronig
        self.gaussian = sigma
        self.amp = A
        self.asymmetry = asym
        self.asymmetryDoniach = asymD
        self.is_singlet = singlet
        self.branching_ratio = branch
        self.spinOrbitSplit = split

        self.amp_2 = self.branching_ratio*self.amp*(self.lorentz + self.gaussian)/(self.lorentz_CK + self.gaussian)
        

        """
        Calculate the Voigt lineshape with variable peak position and intensity using a convolution of a Gaussian and a Lorentzian distribution.

        Parameters:
            x (array-like): The input array of independent variables.
            self.bindingEnergy (float): The position of the peak.
            self.gaussian (float): The standard deviation of the Gaussian distribution.
            self.lorentz (float): The full width at half maximum (FWHM) of the Lorentzian distribution.
            intensity (float): The intensity of the peak.

        Returns:
            array-like: The values of the Voigt lineshape at the given x-values.
        """

        '''
        
        x_array = np.arange(269,295,0.05) #Stepsize of 0.05 eV

        y_array_voigt = np.zeros(len(x_array))

        # creating some noise to add the the y-axis data
        #y_noise_voigt = np.random.normal(0, 0.6, len(x)) #Random noise with standard divation = 0.6 which is about the gaussian ditribution added from our XPS system
        #y_array_voigt += y_noise_voigt


        def fake_data(bindingEnergy, sigma, lorentz, amp):
            lorentz_CK = 0.5625
            branching_ratio = 0.5
            br_CK = 0.5
            spinOrbitSplit = -6.5
            asym = 5.0
            data_range= max(x_array) - min(x_array)
            data_range /= 2
            middle = min(x_array)+data_range
            num_points = len(x_array)
            offset = bindingEnergy-middle
            new_sigma = sigma/(2*np.sqrt(2*np.log(2)))
            x_values, dx = np.linspace(-data_range,data_range,num_points,retstep=True)

            HWHM = lorentz/2
            lorentzLeft = HWHM*asym
            HWHM_CK = lorentz_CK/2
            lorentzLeft_CK = HWHM_CK*asym
            numP = len(x_array)
            yDoubleL = [0]*numP
            #lorentzian = (lorentz / (2*np.pi*(np.power(x_values+offset, 2) + np.power(lorentz/2, 2))))
            
            for i in range(len(x)):
                        
                if round(x_values[i] - offset, 5) <= 0:

                    yDoubleL[i] = 1 / ( 1 + np.power( (x_values[i] - offset)/HWHM, 2 ) ) / np.pi
                else:

                    yDoubleL[i] = 1 / ( 1 + np.power( (x_values[i] - offset)/lorentzLeft, 2 ) ) / np.pi
            lorentzian = yDoubleL

            """
            nu1 = lorentz
            de11 = np.power((x_values-offset), 2) + np.power((lorentz/2), 2)
            de12 = 2*np.pi*de11
            lorentzian1 = (nu1/de12)/(1-br_CK)

            nu2 = br_CK*lorentz_CK
            de21 = np.power(x_values-offset+spinOrbitSplit,2) + np.power(lorentz_CK/2,2)
            de22 = 2*np.pi*de21
            lorentzian2 = (nu2/de22)/(1-br_CK)
            lorentzian = lorentzian1 + lorentzian2
            

            """
            gaussian = np.exp(-np.power(x_values, 2) / (2*(np.power(new_sigma, 2)))) / (new_sigma * np.sqrt(2 * np.pi))
            
            

            final_voigt = scipy.signal.fftconvolve(gaussian,lorentzian,'same')
            #final_voigt = final_voigt / (HWHM/2 + lorentzLeft/2)
            #final_voigt = final_voigt / (HWHM/4 + lorentzLeft/4 + HWHM_CK/4 + lorentzLeft_CK/4)
        
            scale = max(final_voigt)
            for i in range(len(final_voigt)):
                final_voigt[i] *= (amp/scale)
            return final_voigt

        peak_1 = fake_data(279.5, 1.0, 0.25, 15000) #adventitious carbon peak energy
        peak_2 = fake_data(286, 1, 0.5625, 4957.99058) #290 = 280 idk why it shifts it this way (higher is lower and lower is higher)
        #peak_3 = fake_data(280, 1.2, 0.2, 375)
        y_array_voigt += peak_1
        y_array_voigt += peak_2
        #y_array_voigt += peak_3
        #std = 0.1 * np.std(y_array_voigt) #10% gaussian noise --> Typical of our XPS system
        std = 0
        y_noise_voigt = np.random.normal(1000, std, len(x_array))
        y_array_voigt += y_noise_voigt
        x_array = x_array[::-1]
        y_array_voigt = y_array_voigt[::-1]
        print("HELLO")
        for i in range(len(x_array)):
            print(x_array[i], y_array_voigt[i])
        '''
        
        
        
        
  
        
        
        
        

        
        
        

        if self.gaussian ==0:
            if self.paramRange == '': #container mode has no ranges (Phase 3d)
                self.gaussian += 0.01
            else:
                self.gaussian += self.paramRange['Gaussian'][2]

        self.sigma = self.gaussian/(2*np.sqrt(2*np.log(2)))

        data_range= max(x) - min(x)
        data_range /= 2
        middle = min(x)+data_range
        offset = self.bindingEnergy-middle
        

        num_points = len(x)
        x_values, dx = np.linspace(-data_range,data_range,num_points,retstep=True)

        #lorentz
        #                        gamma
        #            ----------------------------
        #                           2    gamma   2
        #            2 Pi ((x - BE ) + ( -----  )  )
        #                                  2
        lorentzian = (self.lorentz / (2*np.pi*(np.power(x_values+offset, 2) + np.power(self.lorentz/2, 2))))

        #gaussian
        #                    2
        #        - (x - BE)
        #        ------------
        #                2
        #         2 sigma
        #    e
        #    -----------------
        #            _____
        #    sigma  |/2 Pi

        #gaussian = np.exp(-np.power(x_values, 2) / (2*(np.power(self.gaussian, 2)))) / (self.gaussian * np.sqrt(2 * np.pi))
        gaussian = np.exp(-np.power(x_values, 2) / (2*(np.power(self.sigma, 2)))) / (self.sigma * np.sqrt(2 * np.pi))

        #Alternative convolution method: Slow even with nb.jit
        '''
        lorentz = lorentzian
        gauss = gaussian
        L, G = len(lorentz), len(gauss)
        cen = G // 2
        conv = np.zeros(L)
        for l in range(L):
            i_min = max(0, cen - l)
            i_max = min(G, L + cen - l)
            temp_conv = np.dot(gauss[i_min:i_max], lorentz[l - cen + i_min : l - cen + i_max])
            conv[l] = temp_conv
        final_voigt = conv
        '''
        final_voigt = scipy.signal.fftconvolve(gaussian,lorentzian,'same')

        if self.is_singlet == False: #If it is a doublet

            #KE Condition
            if x[0] < x[-1]:
                if self.is_coster_kronig == True: #Should symmetrical(not double lorentzian) peaks be allowed to be coster-kronig? 
                    nu1 = self.lorentz
                    de11 = np.power((x_values+offset), 2) + np.power((self.lorentz/2), 2)
                    de12 = 2*np.pi*de11
                    lorentzian1 = (nu1/de12)#/(1-self.branching_ratio)

                    nu2 = self.lorentz_CK#*self.branching_ratio
                    de21 = np.power(x_values+offset-self.spinOrbitSplit,2) + np.power(self.lorentz_CK/2,2)
                    de22 = 2*np.pi*de21
                    lorentzian2 = (nu2/de22)#/(1-self.branching_ratio)
                    lorentzian = lorentzian1 #+ lorentzian2  
                    


                else:

                    nu1 = self.lorentz
                    de11 = np.power((x_values+offset), 2) + np.power((self.lorentz/2), 2)
                    de12 = 2*np.pi*de11
                    lorentzian1 = (nu1/de12)/(1-self.branching_ratio)

                    nu2 = self.branching_ratio*self.lorentz
                    de21 = np.power(x_values+offset-self.spinOrbitSplit,2) + np.power(self.lorentz/2,2)
                    de22 = 2*np.pi*de21
                    lorentzian2 = (nu2/de22)/(1-self.branching_ratio)
                    lorentzian = lorentzian1 + lorentzian2
            
            else: #BE

                if self.is_coster_kronig == True:

                    nu1 = self.lorentz
                    de11 = np.power((x_values+offset), 2) + np.power((self.lorentz/2), 2)
                    de12 = 2*np.pi*de11
                    lorentzian1 = (nu1/de12)#/(1-self.branching_ratio)

                    nu2 = self.lorentz_CK#*self.branching_ratio
                    de21 = np.power(x_values+offset+self.spinOrbitSplit,2) + np.power(self.lorentz_CK/2,2)
                    de22 = 2*np.pi*de21
                    lorentzian2 = (nu2/de22)#/(1-self.branching_ratio)
                    lorentzian = lorentzian1 #+ lorentzian2

               
                    final_voigt = scipy.signal.fftconvolve(gaussian,lorentzian,'same')
                else:

                    nu1 = self.lorentz
                    de11 = np.power((x_values+offset), 2) + np.power((self.lorentz/2), 2)
                    de12 = 2*np.pi*de11
                    lorentzian1 = (nu1/de12)/(1-self.branching_ratio)

                    nu2 = self.branching_ratio*self.lorentz
                    de21 = np.power(x_values+offset+self.spinOrbitSplit,2) + np.power(self.lorentz/2,2)
                    de22 = 2*np.pi*de21
                    lorentzian2 = (nu2/de22)/(1-self.branching_ratio)
                    lorentzian = lorentzian1 + lorentzian2

            final_voigt = scipy.signal.fftconvolve(gaussian,lorentzian,'same')


        #normalize the height so that intensity is the height of the max of the peak
        
        scale = max(final_voigt)
        for i in range(len(final_voigt)):
            final_voigt[i] *= (self.amp/scale)
        if x[0] < x[-1]:
            final_voigt = final_voigt[::-1]
        else:
            pass

        #Combine two singlets for Coster-Kronig
        if self.is_coster_kronig == True:
            if lorentzian2[0] == 0:
                pass 
            else:
                voigt_CK = scipy.signal.fftconvolve(gaussian,lorentzian2,'same')
                #doubleLorentz_CK = doubleLorentz_CK / (HWHM_CK/2 + lorentzLeft_CK/2)
                scale2 = max(voigt_CK)
                for i in range(len(voigt_CK)):
                    voigt_CK[i] *= (self.amp_2/scale2) #Using 2nd peak amp --> Find this, then relate to lorentz_CK

                #Combine for final y values
                final_voigt = final_voigt + voigt_CK
                
        
        #returns, but also updates the yValues of the fit to improve efficiency, we can call that instead of recalculating every time
        self.peak_y = final_voigt
        #peak.voigt = voigt

        return final_voigt





    def gaussFunc(self,x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split):


        #Check to see if the data is in KE energy
        self.bindingEnergy = BE
        self.lorentz = width
        self.lorentz_CK = width_CK
        self.is_coster_kronig = coster_kronig
        self.gaussian = sigma
        self.amp = A
        self.asymmetry = asym
        self.asymmetryDoniach = asymD
        self.is_singlet = singlet
        self.branching_ratio = branch
        self.spinOrbitSplit = split
     

        if self.gaussian ==0:
            if self.paramRange == '': #container mode has no ranges (Phase 3d)
                self.gaussian += 0.01
            else:
                self.gaussian += self.paramRange['Gaussian'][2]

        self.sigma = self.gaussian/(2*np.sqrt(2*np.log(2)))


        data_range= max(x) - min(x)
        data_range /= 2

        middle = min(x)+data_range
        offset = self.bindingEnergy-middle
        num_points = len(x)
        x_values, dx = np.linspace(-data_range,data_range,num_points,retstep=True)



        #gaussian = np.exp(-np.power(x_values+offset, 2) / (2*(np.power(self.gaussian, 2)))) / (self.gaussian * np.sqrt(2 * np.pi))
        gaussian = np.exp(-np.power(x_values+offset, 2) / (2*(np.power(self.sigma, 2)))) / (self.sigma * np.sqrt(2 * np.pi))

        #gaussian = np.exp(-np.power(x - self.bindingEnergy, 2) / (2 * np.power(self.gaussian, 2))) / (self.gaussian * np.sqrt(2 * np.pi))

        scale = max(gaussian)
        for i in range(len(gaussian)):
            gaussian[i] *= (self.amp/scale)

        #self.peak_y = gaussian

        if self.is_singlet == False: #If it is a doublet
            if x[0] < x[-1]:
                #gaussian1 = (np.exp(-np.power(x_values+offset, 2) / (2 * np.power(self.gaussian, 2))) / (self.gaussian * np.sqrt(2 * np.pi)))/(1-self.branching_ratio)
                #gaussian2 = self.branching_ratio * (np.exp(-np.power((x_values + (offset - self.spinOrbitSplit)), 2) / (2 * np.power(self.gaussian, 2))) / (self.gaussian * np.sqrt(2 * np.pi)))/(1-self.branching_ratio)
                gaussian1 = (np.exp(-np.power(x_values+offset, 2) / (2 * np.power(self.sigma, 2))) / (self.sigma * np.sqrt(2 * np.pi)))/(1-self.branching_ratio)
                gaussian2 = self.branching_ratio * (np.exp(-np.power((x_values + (offset - self.spinOrbitSplit)), 2) / (2 * np.power(self.sigma, 2))) / (self.sigma * np.sqrt(2 * np.pi)))/(1-self.branching_ratio)
            else:
                #gaussian1 = (np.exp(-np.power(x_values+offset, 2) / (2 * np.power(self.gaussian, 2))) / (self.gaussian * np.sqrt(2 * np.pi)))/(1-self.branching_ratio)
                #gaussian2 = self.branching_ratio * (np.exp(-np.power((x_values + (offset + self.spinOrbitSplit)), 2) / (2 * np.power(self.gaussian, 2))) / (self.gaussian * np.sqrt(2 * np.pi)))/(1-self.branching_ratio)
                gaussian1 = (np.exp(-np.power(x_values+offset, 2) / (2 * np.power(self.sigma, 2))) / (self.sigma * np.sqrt(2 * np.pi)))/(1-self.branching_ratio)
                gaussian2 = self.branching_ratio * (np.exp(-np.power((x_values + (offset + self.spinOrbitSplit)), 2) / (2 * np.power(self.sigma, 2))) / (self.sigma * np.sqrt(2 * np.pi)))/(1-self.branching_ratio)
            gaussian = gaussian1 + gaussian2


            scale = max(gaussian)
            for i in range(len(gaussian)):
                gaussian[i] *= (self.amp/scale)

        '''
        if x[0] < x[-1]: #KE CONDITION
            gaussian = gaussian[::-1]
        else:
            pass
        '''
        self.peak_y = gaussian
        

        return gaussian





    def lorentzFunc(self,x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split):

        global bindingEnergy
        self.bindingEnergy = BE
        self.lorentz = width
        self.lorentz_CK = width_CK
        self.is_coster_kronig = coster_kronig
        self.gaussian = sigma
        self.amp = A
        self.asymmetry = asym
        self.asymmetryDoniach = asymD
        self.is_singlet = singlet
        self.branching_ratio = branch
        self.spinOrbitSplit = split



        data_range= max(x) - min(x)
        data_range /= 2

        middle = min(x)+data_range
        offset = self.bindingEnergy-middle
        num_points = len(x)
        x_values, dx = np.linspace(-data_range,data_range,num_points,retstep=True)
        #lorentz
        

       
        #self.peak_y = lorentzian
        if self.is_singlet == True:
            lorentzian = (self.lorentz / (2*np.pi*(np.power(x_values+offset, 2) + np.power(self.lorentz/2, 2))))
            scale = max(lorentzian)
            for i in range(len(lorentzian)):
                lorentzian[i] *= (self.amp/scale)




        elif self.is_singlet == False: #If it is a doublet
            #No clue why but self.lorentz does not need to be made into the HWHM for the doublet. Branching ratio is still off for some reason.
            if x[0] < x[-1]: #KE 

                    nu1 = self.lorentz
                    de11 = np.power((x_values+offset), 2) + np.power((self.lorentz/2), 2)
                    de12 = 2*np.pi*de11
                    lorentzian1 = (nu1/de12)/(1-self.branching_ratio)

                    nu2 = self.branching_ratio*self.lorentz
                    de21 = np.power(x_values+offset-self.spinOrbitSplit,2) + np.power(self.lorentz/2,2)
                    de22 = 2*np.pi*de21
                    lorentzian2 = (nu2/de22)/(1-self.branching_ratio)

                    lorentzian = (lorentzian1 + lorentzian2)


                    scale = max(lorentzian)
                    for i in range(len(lorentzian)):
                        lorentzian[i] *= (self.amp/scale)
            else: #BE

                    nu1 = self.lorentz
                    de11 = np.power((x_values+offset), 2) + np.power((self.lorentz/2), 2)
                    de12 = 2*np.pi*de11
                    lorentzian1 = (nu1/de12)/(1-self.branching_ratio)

                    nu2 = self.branching_ratio*self.lorentz
                    de21 = np.power(x_values+offset+self.spinOrbitSplit,2) + np.power(self.lorentz/2,2)
                    de22 = 2*np.pi*de21
                    lorentzian2 = (nu2/de22)/(1-self.branching_ratio)

                    lorentzian = (lorentzian1 + lorentzian2)


                    scale = max(lorentzian)
                    for i in range(len(lorentzian)):
                        lorentzian[i] *= (self.amp/scale)
        
        '''
        if x[0] < x[-1]: #KE CONDITION
            lorentzian = lorentzian[::-1]
        else:
            pass
        '''

        self.peak_y = lorentzian


        return lorentzian






    def doubleLorentzFunc(self,x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split):
        #Same as voigt but with an added asymmetry factor to the width (lorentz)in the lorentzian equation
        global bindingEnergy
        self.bindingEnergy = BE
        self.lorentz = width
        self.lorentz_CK = width_CK
        self.asym_CK = asym_CK
        self.is_coster_kronig = coster_kronig
        self.gaussian = sigma
        self.amp = A
        self.asymmetry = asym
        self.asymmetryDoniach = asymD
        self.is_singlet = singlet
        self.branching_ratio = branch
        self.spinOrbitSplit = split

        #self.amp_2 = 0.4*self.amp #Amp of second peak 
        LFWHM = (self.lorentz*self.asymmetry)/2 + self.lorentz/2
        LFWHM_CK = (self.lorentz_CK*self.asym_CK)/2 + self.lorentz_CK/2 #Used to find correct amp value that has the right area ratio
        #LFWHM_CK = (self.lorentz_CK*self.asymmetry)/2 + self.lorentz_CK/2

        self.amp_2 = self.branching_ratio*self.amp*(LFWHM + self.gaussian)/(LFWHM_CK + self.gaussian) #estimate needs the LFWHM which is different on each side of peak
        #self.amp_2 = self.branching_ratio*self.amp*(LFWHM*self.asymmetry + self.gaussian)/(LFWHM_CK*self.asym_CK + self.gaussian)
        
        # Calculate the Gaussian component
        if self.gaussian == 0:
            self.gaussian = .01

        self.sigma = self.gaussian/(2*np.sqrt(2*np.log(2)))
        


          
        data_range= max(x) - min(x)
        data_range /= 2

        middle = min(x)+data_range
        #offset = self.bindingEnergy-middle
    
        #sos_offset = offset - self.spinOrbitSplit//2 #offset for spin-orbit splitting
        #sos_offset = -sos_offset
        #print(offset, sos_offset)
        
        num_points = len(x)
        x_values, dx = np.linspace(-data_range,data_range,num_points,retstep=True)
        dx = round(dx, 5)

        #common multipliers
     
        self.bindingEnergy = self.bindingEnergy #+ dx*mult
        offset = self.bindingEnergy-middle #- offset_adj#- dx*mult #Data is offset by dx*val for some reason --> DL property. Mentioned by Alberto.

        #gaussian = np.exp(-np.power(x - self.bindingEnergy, 2) / (2 * np.power(self.gaussian, 2))) / (self.gaussian * np.sqrt(2 * np.pi))
        #gaussian = np.exp(-np.power(x_values, 2) / ((np.power(self.gaussian, 2)))) / (self.gaussian * np.sqrt(np.pi)) #Took away 2*np.power(self.gaussian,2) and np.sqrt(2*np.pi)
        #gaussian = np.exp(-np.power(x_values, 2) / ((np.power(self.sigma, 2)))) / (self.sigma * np.sqrt(np.pi))
        gaussian = np.exp(-np.power(x_values, 2) / (2*(np.power(self.sigma, 2)))) / (self.sigma * np.sqrt(2 * np.pi))
        # Calculate the Lorentzian component

        #Added +offset to all Lorentzian functions instead of Gauss --> It has to be inside the equation to make the leftside of the peak more lorentzian
        new_x = np.zeros(len(x))
        numP = len(x)
        HWHM = self.lorentz/2
        lorentzLeft = HWHM*self.asymmetry #Width of left side of peak due to asymmetry

        HWHM_CK = self.lorentz_CK/2
        lorentzLeft_CK = HWHM_CK*self.asymmetry
        
        lorentzian2 = 0

        yDoubleL = [0]*numP
        yDoubleL_CK = [0]*numP

        if self.is_singlet == True:
            for i in range(len(x)):

                #KE Condition
                if x[0] < x[-1]:
                    if round(x_values[i] - offset, 5) <= 0:

                        yDoubleL[i] = 1 / ( 1 + np.power( (x_values[i] - offset)/lorentzLeft, 2 ) ) / np.pi
                    else:

                        yDoubleL[i] = 1 / ( 1 + np.power( (x_values[i] - offset)/HWHM, 2 ) ) / np.pi


                else:
                    if round(x_values[i] + offset,5) <= 0:

                        yDoubleL[i] = 1 / ( 1 + np.power( (x_values[i] + offset)/lorentzLeft, 2 ) ) / np.pi
                    else:

                        yDoubleL[i] = 1 / ( 1 + np.power( (x_values[i] + offset)/HWHM, 2 ) ) / np.pi

                #yDoubleL[i] = yDoubleL[i] / (lorentzLeft/2 + HWHM/2)

            lorentzian = yDoubleL

        elif self.is_singlet == False: #If it is a doublet
            #KE Condition
            if x[0] < x[-1]: #Why is the offset different if it is a singlet or doublet for KE?
                if self.is_coster_kronig == True: #If it has coster-kronig effects --> This applies to all peaks. Hot to make it see individual peaks?
                    for i in range(len(x)):

                        #Low BE peak
                        if round(x_values[i] - offset, 5) <= 0:

                            yDoubleL[i] = 1 / ( 1 + np.power( (x_values[i] - offset)/lorentzLeft, 2 ) ) / np.pi
                        else:

                            yDoubleL[i] = 1 / ( 1 + np.power( (x_values[i] - offset)/HWHM, 2 ) ) / np.pi

                        #High BE peak
                        if round(x_values[i] - offset + self.spinOrbitSplit, 5) <= 0:

                            yDoubleL_CK[i] = 1 / ( 1 + np.power( (x_values[i] + (-offset + self.spinOrbitSplit))/lorentzLeft_CK, 2 ) ) / np.pi
                        else:

                            yDoubleL_CK[i] = 1 / ( 1 + np.power( (x_values[i] + (-offset + self.spinOrbitSplit))/HWHM_CK, 2 ) ) / np.pi

                    lorentzian = yDoubleL
                    lorentzian2 = yDoubleL_CK



                else: #If not coster-kronig (just doublet) 

                    for i in range(len(x)):

                        if round(x_values[i] - offset,5) <= 0:
                            yDoubleL_right = 1 / (1 + np.power( (x_values[i] - offset)/ lorentzLeft, 2)) #/ (1 - self.branching_ratio) / np.pi / lorentzLeft/2

                        else:
                            
                            yDoubleL_right = 1 / (1 + np.power( (x_values[i] - offset)/ HWHM, 2)) #/ (1 - self.branching_ratio) / np.pi / HWHM/2


                        #Higher BE peak
                        if round(x_values[i] - offset + self.spinOrbitSplit, 5) <= 0:
                            yDoubleL_left = self.branching_ratio / (1 + np.power( (x_values[i] + (-offset + self.spinOrbitSplit)) / lorentzLeft, 2)) #/ (1 - self.branching_ratio) / np.pi / lorentzLeft/2
                        else:
                            
                            yDoubleL_left = self.branching_ratio / (1 + np.power( (x_values[i] + (-offset + self.spinOrbitSplit)) / HWHM, 2)) #/ (1 - self.branching_ratio) / np.pi / HWHM/2

                        yDoubleL[i] = (yDoubleL_right + yDoubleL_left) / (1 - self.branching_ratio) / np.pi



            #BE
            else:
                if self.is_coster_kronig == True: #If it has coster-kronig effects --> This applies to all peaks. Hot to make it see individual peaks?
                    for i in range(len(x)):

                        #Low BE peak
                        if round(x_values[i] + offset, 5) <= 0:

                            yDoubleL[i] = 1 / ( 1 + np.power( (x_values[i] + offset)/lorentzLeft, 2 ) ) / np.pi
                        else:

                            yDoubleL[i] = 1 / ( 1 + np.power( (x_values[i] + offset)/HWHM, 2 ) ) / np.pi

                        #High BE peak
                        if round(x_values[i] + offset + self.spinOrbitSplit, 5) <= 0:

                            yDoubleL_CK[i] = 1 / ( 1 + np.power( (x_values[i] + (offset + self.spinOrbitSplit))/lorentzLeft_CK, 2 ) ) / np.pi
                        else:

                            yDoubleL_CK[i] = 1 / ( 1 + np.power( (x_values[i] + (offset + self.spinOrbitSplit))/HWHM_CK, 2 ) ) / np.pi

                    lorentzian = yDoubleL
                    lorentzian2 = yDoubleL_CK


                else: #If not coster-kronig (just doublet)

                    for i in range(len(x)):

                        if round(x_values[i] + offset,5) <= 0:
                            yDoubleL_right = 1 / (1 + np.power( (x_values[i] + offset)/ lorentzLeft, 2)) #/ (1 - self.branching_ratio) / np.pi / lorentzLeft/2

                        else:
                            
                            yDoubleL_right = 1 / (1 + np.power( (x_values[i] + offset)/ HWHM, 2)) #/ (1 - self.branching_ratio) / np.pi / HWHM/2


                        #Higher BE peak
                        if round(x_values[i] + offset + self.spinOrbitSplit, 5) <= 0:
                            yDoubleL_left = self.branching_ratio / (1 + np.power( (x_values[i] + (offset + self.spinOrbitSplit)) / lorentzLeft, 2)) #/ (1 - self.branching_ratio) / np.pi / lorentzLeft/2
                        else:
                            
                            yDoubleL_left = self.branching_ratio / (1 + np.power( (x_values[i] + (offset + self.spinOrbitSplit)) / HWHM, 2)) #/ (1 - self.branching_ratio) / np.pi / HWHM/2

                        yDoubleL[i] = (yDoubleL_right + yDoubleL_left) / (1 - self.branching_ratio) / np.pi

            lorentzian = yDoubleL

        doubleLorentz = scipy.signal.fftconvolve(gaussian,lorentzian,'same')

      
        #if self.is_singlet == False:
            #doubleLorentz = doubleLorentz / (HWHM/4 + lorentzLeft/4 + HWHM_CK/4 + lorentzLeft_CK/4)
            #doubleLorentz = doubleLorentz / (HWHM/2 + lorentzLeft/2)
        #else:
        #    doubleLorentz = doubleLorentz / (HWHM/2 + lorentzLeft/2)
            
        scale = max(doubleLorentz)
        for i in range(len(doubleLorentz)):
            doubleLorentz[i] *= (self.amp/scale) #got rid of scale

        #Combine two singlets
        if self.is_coster_kronig == True:
            if lorentzian2 == 0:
                pass 
            else:
                doubleLorentz_CK = scipy.signal.fftconvolve(gaussian,lorentzian2,'same')
                #doubleLorentz_CK = doubleLorentz_CK / (HWHM_CK/2 + lorentzLeft_CK/2)
                scale2 = max(doubleLorentz_CK)
                for i in range(len(doubleLorentz_CK)):
                    doubleLorentz_CK[i] *= (self.amp_2/scale2) #Using 2nd peak amp --> Find this, then relate to lorentz_CK

                #Combine for final y values
                doubleLorentz = doubleLorentz + doubleLorentz_CK
        


        '''
        if BE_entry == True: #KE CONDITION
            #doubleLorentz = doubleLorentz[::-1]
            x = x[::-1]
        '''
       

        #print("HELLO", 2*HWHM, 2*HWHM_CK)
        self.peak_y = doubleLorentz
        #for i in np.arange(1, numP, 1):
            #print(x[i], doubleLorentz[i])
       
        return doubleLorentz





    def doniachSunjicFunc(self,x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split):
        #Formula for Doniach-Sunjic peak equation:
        #
        #             cos(( pi  *  alpha ) / 2  +  (1 -  alpha ) * arctan((x - center) /  self.lorentz ))
        # Func =   ----------------------------------------------------------------------------
        #                                2          2  (1 - alpha) / 2
        #                  ( (x - center)   +   gamma )
        global bindingEnergy
        self.bindingEnergy = BE
        self.lorentz = width
        self.lorentz_CK = width_CK
        self.is_coster_kronig = coster_kronig
        self.gaussian = sigma
        self.amp = A
        self.asymmetry = asym
        self.asymmetryDoniach = asymD
        self.is_singlet = singlet
        self.branching_ratio = branch
        self.spinOrbitSplit = split
    

        data_range= max(x) - min(x)
        data_range /= 2

        middle = min(x)+data_range
        offset = self.bindingEnergy-middle
        num_points = len(x)
        x_values, dx = np.linspace(-data_range,data_range,num_points,retstep=True)

        if self.gaussian ==0:
            if self.paramRange == '': #container mode has no ranges (Phase 3d)
                self.gaussian += 0.01
            else:
                self.gaussian += self.paramRange['Gaussian'][2]

        self.sigma = self.gaussian/(2*np.sqrt(2*np.log(2)))

        #gaussian = np.exp(-np.power(x_values, 2) / ((np.power(self.gaussian, 2)))) / (self.gaussian * np.sqrt(np.pi)) #Took away 2*np.power(self.gaussian,2) and np.sqrt(2*np.pi)
        #gaussian = np.exp(-np.power(x_values, 2) / ((np.power(self.sigma, 2)))) / (self.sigma * np.sqrt(np.pi))
        gaussian = np.exp(-np.power(x_values, 2) / (2*(np.power(self.sigma, 2)))) / (self.sigma * np.sqrt(2 * np.pi))
        HWHM = self.lorentz/2
        if self.is_singlet == True:
            #self.asymmetryDoniach = -self.asymmetryDoniach #Opposite sign because it is in BE not KE


            #Singlet Doniach-Sunjic Equation:
            cos1 = (np.pi*self.asymmetryDoniach)/2
            cos2 = (1-self.asymmetryDoniach)*np.arctan((x_values+offset)/HWHM)
            numerator = np.cos(cos1 + cos2)
            de1 = pow(x_values+offset, 2) + pow(HWHM, 2)
            powDe = (self.asymmetryDoniach-1)/2
            denominator = pow(de1, powDe)
            lorentzian = numerator*denominator



        #Doesnt give doublet right now. Do not know why --> it follows the same methods as Aanalyzer
        elif self.is_singlet == False: #If it is a doublet

                cos1 = (np.pi*self.asymmetryDoniach)/2
                cos2 = (1-self.asymmetryDoniach)*np.arctan((x_values+offset)/HWHM)
                cos3 = (1-self.asymmetryDoniach)*np.arctan((x_values+self.spinOrbitSplit + offset)/HWHM) #doublet only for the second part of the equation
                numerator1 = np.cos(cos1 + cos2)
                numerator2 = np.cos(cos1 + cos3)
                de1 = pow(x_values+offset, 2) + pow(HWHM, 2)
                de2 = pow(x_values + offset + self.spinOrbitSplit, 2) + pow(HWHM, 2)
                powDe = (self.asymmetryDoniach-1)/2
                denominator1 = pow(de1, powDe)
                denominator2 = pow(de2, powDe)
                lorentzian = numerator1*denominator1 + self.branching_ratio*(numerator2*denominator2)




        doniachSunjic = scipy.signal.fftconvolve(gaussian,lorentzian,'same')

        #normalize the height so that intensity is the height of the max of the peak
        scale = max(doniachSunjic)
        for i in range(len(doniachSunjic)):
            doniachSunjic[i] *= (self.amp/scale)


        if x[0] < x[-1]: #KE CONDITION
            doniachSunjic = doniachSunjic[::-1]
        else:
            pass

        self.peak_y = doniachSunjic

        return doniachSunjic





#---------------------------Backgrounds Class---------------------------------------------------------------


class background(peak):
    """
    Background class is largely the same format as the peaks class,getY, set,get, function largely the same without
    most of the hassle the peaks come with

    One thing of note, most backgrounds are held by the individual, but Shirley backgrounds are held by a peak associated with them,
    and the peak takes on the management of it,but inside the background class it is the same as any other background

    To add a Background type:
    Add its bkgnType to the picker in the GUI
    Give its parameters ranges and pass them into the write_ini in the gui
    add the range to ini_parser and also the paramrange dict in xps_neo
    add that bkgnType in init here, and assign the function to self.bkgn
    add the get function for it
    add the set function in the gui xps_fit
    add the bgknType option to the self.bkgn_types at the top of xps_analysis2

    """
    def __init__(self,paramRange,bkgnType, peakType):
        self.bkgnType = bkgnType
        self.paramRange= paramRange
        self.peakType = peakType

        if paramRange == '':
            # Analysis container mode (Phase 3d): values are assigned via
            # the set_* methods afterwards; no ranges, no RNG draws.
            self.bindingEnergy = 0.0
            self.k = 0.0
            self.background = 0.0
            self.slope = 0.0
            self.baseline_value = 0.0
            if self.bkgnType == 'SVSC':
                self.bkgn = self.shirley_bkgn_again
            elif self.bkgnType == 'Shirley':
                self.bkgn = self.shirley
            elif self.bkgnType.lower() == 'linear':
                self.bkgn = self.linear_background
            elif self.bkgnType == 'Exponential':
                self.bkgn = self.exponential_bkgn
            elif self.bkgnType == 'Polynomial 1':
                self.bkgn = self.polynomial1
            elif self.bkgnType == 'Polynomial 2':
                self.bkgn = self.polynomial2
            elif self.bkgnType == 'Polynomial 3':
                self.bkgn = self.polynomial3
            elif self.bkgnType == '3-Param Tougaard':
                self.CTou3 = 1000.0 #Initial inputs in Aanalyzer
                self.DTou3 = 13300.0
                self.bkgn = self.Tougaard3Param
            elif self.bkgnType == '2-Param Tougaard':
                self.bkgn = self.Tougaard2Param
            elif self.bkgnType == 'Baseline':
                self.bkgn = self.baseline2
            else:
                print("Error Choosing Background in init of xps_fit")
                print("Background read as: " + str(self.bkgnType))
                sys.exit(1)
            self.yBkgn = []
            return

        peak_class = peak(self.paramRange, self.peakType)
        BE_input = np.min(self.paramRange['BE'])
        BE_adjustment = peak_class.bindingEnergy
        self.bindingEnergy = BE_input + BE_adjustment
        








        if self.bkgnType == 'SVSC':
            #self.bkgn = self.shirley_Sherwood

            #New Shirley Background method that uses BE location of each peak along with the HWHM (lorentz/2)
           
            self.bkgn = self.shirley_bkgn_again
           

            self.k_range = np.arange(paramRange['k_range'][0],paramRange['k_range'][1],paramRange['k_range'][2])
            self.k = np.random.choice(self.k_range)
            
            
           
            #self.backgroundShirRange = np.arange(paramRange['Shirley Background'][0],paramRange['Shirley Background'][1],paramRange['Shirley Background'][2])
            #self.backgroundShirley = np.random.choice(self.backgroundShirRange)
        elif self.bkgnType == 'Shirley':
            #self.bkgn = self.shirley_Sherwood
            
            self.bkgn = self.shirley
            

            self.k_range = np.arange(paramRange['k_range'][0],paramRange['k_range'][1],paramRange['k_range'][2])
            self.k = np.random.choice(self.k_range)

        elif self.bkgnType.lower() == 'linear':
            self.bkgn = self.linear_background
            self.backgroundRange = np.arange(paramRange['Background'][0],paramRange['Background'][1],paramRange['Background'][2])
            self.slopeRange = np.arange(paramRange['Slope'][0],paramRange['Slope'][1],paramRange['Slope'][2])

            #self.background is the b value in y = mx+b
            self.background = np.random.choice(self.backgroundRange)

            self.slope = np.random.choice(self.slopeRange)
            #self.slope = 0
            #
        elif self.bkgnType == 'Exponential':
            #self.bkgn = self.exponential_bkgn
            self.bkgn = self.exponential_bkgn
            #self.bkgn = self.new_exponential

        elif self.bkgnType == 'Polynomial 1':
            self.bkgn = self.polynomial1
        elif self.bkgnType == 'Polynomial 2':
            self.bkgn = self.polynomial2
        elif self.bkgnType == 'Polynomial 3':
            self.bkgn = self.polynomial3
        elif self.bkgnType == '3-Param Tougaard':
            self.CTou3 = 1000.0 #Initial inputs in Aanalyzer
            self.DTou3 = 13300.0
            self.bkgn = self.Tougaard3Param
        elif self.bkgnType == '2-Param Tougaard':
            self.bkgn = self.Tougaard2Param
                
        elif self.bkgnType == 'Baseline':
            self.bkgn = self.baseline2

            self.baselineRange = np.arange(paramRange['baseline'][0],paramRange['baseline'][1],paramRange['baseline'][2])
            self.baseline_value = np.random.choice(self.baselineRange)
        else:
            print("Error Choosing Background in init of xps_fit")
            print("Background read as: " + str(self.bkgnType))
            sys.exit(1)
        self.yBkgn = []

    _PICKLE_SKIP = ("bkgn", "k_range", "backgroundRange", "slopeRange",
                    "baselineRange")

    def __getstate__(self):
        return {k: v for k, v in self.__dict__.items()
                if k not in self._PICKLE_SKIP}

    def __setstate__(self, state):
        self.__dict__.update(state)
        if self.bkgnType == 'SVSC':
            self.bkgn = self.shirley_bkgn_again
        elif self.bkgnType == 'Shirley':
            self.bkgn = self.shirley
        elif self.bkgnType.lower() == 'linear':
            self.bkgn = self.linear_background
        elif self.bkgnType == 'Exponential':
            self.bkgn = self.exponential_bkgn
        elif self.bkgnType == 'Polynomial 1':
            self.bkgn = self.polynomial1
        elif self.bkgnType == 'Polynomial 2':
            self.bkgn = self.polynomial2
        elif self.bkgnType == 'Polynomial 3':
            self.bkgn = self.polynomial3
        elif self.bkgnType == '3-Param Tougaard':
            self.bkgn = self.Tougaard3Param
        elif self.bkgnType == '2-Param Tougaard':
            self.bkgn = self.Tougaard2Param
        elif self.bkgnType == 'Baseline':
            self.bkgn = self.baseline2

    def getY(self,x,y, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split,scale_var, tot_area):
        y_val = []
        first = y[0]
        if first < 1:
            scale_val = round(1/first)*100
        else:
            scale_val = first
        
        if scale_var == True:
            
            y_val = y/first #Dividing every element by the first value 
            
            y_val = y*scale_val #Multiply by 1000 to scale it
               
        else:
            y_val = y
        y = y_val
      
        self.get_Background(x,y, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, tot_area)
        
        return self.yBkgn



    def get_Background(self,x,y, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, tot_area):
        
        self.bkgn(x,y, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, tot_area)











    def mutate(self,chance):
        if self.bkgnType == 'SVSC':
            self.mutate_k(chance)
        if self.bkgnType == 'Shirley':
            self.mutate_k(chance)
        if self.bkgnType == 'Baseline':
            self.mutate_baseline_value(chance)
        if self.bkgnType == 'Linear': #Added this in for new background --> Got rid of old slope
            self.mutate_background(chance)
            self.mutate_slope(chance)
        '''
        if self.bkgnType == 'Exponential':
            self.mutate_A(chance)
            self.mutate_tau(chance)
        if self.bkgnType == '3-Param Tougaard':
            self.mutate_C(chance)
            self.mutate_D(chance)
        '''


    def mutate_k(self,chance):
        if random.random()*100 < chance:
            self.k = np.random.choice(self.k_range)

    def mutate_baseline_value(self,chance):
        if random.random()*100 < chance:
            self.baseline_value = np.random.choice(self.baselineRange)

    def mutate_background(self,chance):
        if random.random()*100 < chance:
            self.background = np.random.choice(self.backgroundRange)
    def mutate_slope(self,chance):
        if random.random()*100 < chance:
            self.slope = np.random.choice(self.slopeRange)

  


    #Make sure to add in each background here
    def get(self):
       if self.bkgnType == 'SVSC':
            #return [self.k, self.backgroundShirley, self.bkgnType]
            return [self.k, self.bkgnType] #How to make k values for each peak? 
       elif self.bkgnType == 'Shirley':
            return [self.k, self.bkgnType] 
       elif self.bkgnType == 'Linear':
            return [self.background,self.slope,self.bkgnType]
       elif self.bkgnType == 'Exponential':
            return [self.bkgnType]
       elif self.bkgnType == 'Baseline':
            return [self.baseline_value,self.bkgnType]
       elif self.bkgnType == 'Polynomial 1':
            return [self.bkgnType]
       elif self.bkgnType == 'Polynomial 2':
            return [self.bkgnType]
       elif self.bkgnType == 'Polynomial 3':
            return [self.bkgnType]
       elif self.bkgnType == '3-Param Tougaard':
            return [self.bkgnType]
       elif self.bkgnType == '2-Param Tougaard':
            return [self.bkgnType]
    def getType(self):
        return self.bkgnType

    def set_k(self,newVal):
        self.k = newVal

    #def set_backgroundShir(self,newVal):
        #self.backgroundShirley = newVal

    def set_shirley_sherwood(self,params):
        self.k = params[0]
      
        #self.backgroundShirley = params[1]
    
    def set_shirley(self,params):
        self.k = params[0]

    def set_baseline(self,params):
        self.baseline_value = params[0]


    def set_background(self,newVal):
        self.background = newVal

    def set_slope(self,newVal):
        self.slope = newVal

    def set_linear(self,params):
        self.background = params[0]
        self.slope = params[1]


    '''
    def baseline(self,x,y):
        self.y = y
        self.x = x
        data_baseline = peakutils.baseline(x)

        s= data_baseline
        #bkgn_0 = data_baseline.modpoly(y, poly_order=0)[0]
        #funcs = bkgn_0
        self.yBkgn = funcs
        return self.yBkgn
    '''
    def baseline(self,x,y, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, tot_area):
        self.y = y
        self.x = x

         #Check to see if the data is in KE energy


        poly_0 = np.polyfit(x, y, deg=0)
        funcs = np.polyval(poly_0, x)


        return self.yBkgn




    def polynomial1(self,x,y, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, tot_area):
        self.y = y
        self.x = x


        numP = len(x)
        funcs = [0]*numP
        for i in np.arange(0, numP):
            funcs[i] = pow( (x[i]-x[int(np.rint(numP/2))]), 1);

        self.yBkgn = funcs

        return self.yBkgn




    def polynomial2(self,x,y, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, tot_area):
        self.y = y
        self.x = x

        numP = len(x)
        funcs = [0]*numP
        #for i in np.arange(0, numP):
            #funcs[i] = pow( (x[i]-x[int(np.rint(numP/2))]), 2);
        poly_2 = np.polyfit(x, y, deg=2)
        funcs = np.polyval(poly_2, x)

        self.yBkgn = funcs

        return self.yBkgn




    def polynomial3(self,x,y, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, tot_area):
        self.y = y
        self.x = x

        numP = len(x)
        funcs = [0]*numP
        #for i in np.arange(0, numP):
            #funcs[i] = pow( (x[i]-x[int(np.rint(numP/2))]), 3);
        poly_3 = np.polyfit(x, y, deg=3)
        funcs = np.polyval(poly_3, x)

        self.yBkgn = funcs

        return self.yBkgn


    def exponential_bkgn(self,x,y, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, tot_area):
        
        log_y = np.log(y)
        N = 30
        y_left = log_y[:N]
        y_right = log_y[-N:]
        x_left = x[:N]
        x_right = x[-N:]
        y_left_avg = sum(y_left)/N
        y_right_avg = sum(y_right)/N
        x_left_avg = sum(x_left)/N
        x_right_avg = sum(x_right)/N
        m = (y_left_avg - y_right_avg)/abs(x_left_avg - x_right_avg)
        
        b = log_y[0] - m*x[0] 
        
        funcs = np.exp(b)*np.exp(m*x) #+c?

        self.yBkgn = funcs

        return self.yBkgn
        


    '''
    def exponential_bkgn(self,x,y, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split):
        self.y = y
        self.x = x

        #need to make these parameters outside background functions


        numP = len(y)
        funcs =[0]*len(y)

        #Taken from Aanalyzer code --> not sure why exponent is initially set to 1 or 0
        exponent = 1
        deltaExponent = max(abs(exponent / 100), 0.001)
        exponent += deltaExponent


        #Not sure if the x data needs to be flipped for BE instead of KE --> The exponential should be on the left side of the peak not the right
        for j in range(1, numP): #Cut off before numP so the end point is off. Need to fix this in order to scale down the righthand side of the background to the data
            gar = -exponent * (x[j] - x[numP // 2])

            if gar > 30:
                gar = 30
            elif gar < -30:
                gar = -30

            funcs[j] = -(np.exp(gar)) #Added negative sign to flip exponential to be in the -xy plane instead of +xy plane

        self.yBkgn = funcs

        return self.yBkgn
    '''



    def baseline2(self,x,y, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, tot_area):

        bkgn_vals = [self.baseline_value]*len(y)
        self.yBkgn = bkgn_vals

        return self.yBkgn
        return bkgn_vals




    #Integral slope background works for now but is bad. Left side of data is not scaling properly
    #dont know why but it takes a long time to compute the integral slope...
    def slope_bkgn(self,x,y, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, tot_area):
        energyDelay= 0.5 #This is the initial value used in Aanalyzer --> not sure why??
        self.y = y
        self.x = x

        #need to make these parameters outside background functions
        numPeaksUsed = 1
        maxPeaks = 1
        numberBackgrounds = 1
        ma = maxPeaks + numberBackgrounds
        funcs = [0]*len(y)
        numP = len(y)



        for j in np.arange(numP-2, 0, -1): #error changed numP -1 to numP-2
            funcs[j] = (y[j] - y[numP-1]) * (x[j+1] - x[j]) + funcs[j+1] #changed y[numP] to y[numP -1] because of index error

        for j in range(1, numP):
            jDelay = 0
            x_eDelay = x[j] + energyDelay
            while j + jDelay < len(x) and x_eDelay > x[j + jDelay]:
                jDelay += 1
                if j + jDelay > numP:
                    jDelay -=1
                    break

            if j + jDelay > numP-1: #supposed to be just numP --> error
                funcs[j] = 0
            else:
                funcs[j] = funcs[j + jDelay]

        for j in np.arange(numP-2, 0, -1):#error changed numP -1 to numP-2
            funcs[j] = funcs[j] * (x[j+1] - x[j]) + funcs[j+1]


        self.yBkgn = funcs

        return self.yBkgn


    def new_slope(self,x,y, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, tot_area):

        self.x = x
        self.y = y

        N = 10
        y_left = self.y[:N]
        y_right = self.y[N:]
        x_left = self.x[:N]
        x_right = self.x[N:]

        y_left_avg = sum(y_left)/N
        x_left_avg = sum(x_left)/N
        y_right_avg = sum(y_right)/N
        x_right_avg = sum(x_right)/N

        y_slope = (y_left_avg - y_right_avg)/(x_left_avg - x_right_avg)
        y_intercept = self.y[0] - y_slope*self.x[0]
        funcs = y_slope*self.x + y_intercept

        self.yBkgn = funcs

        return funcs


    def linear_background(self,x,y, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, tot_area):
        self.y = y
        self.x = x
        numP = len(self.x)
        self.yBkgn = [0]*numP
        data_range= max(x) - min(x)
        data_range /= 2
        middle = min(x)+data_range
    


        num_points = len(x)
        x_values, dx = np.linspace(-data_range,data_range,num_points,retstep=True)


        #slope = (self.y[-1] - self.y[0])/(self.x[-1] - self.x[0])
        for i in range(numP):
            self.yBkgn[i] = self.linear(self.slope,x_values[i],self.background)
            #print(self.x[i], self.yBkgn[i])
        #print("--------------------------------------------")
        #print(self.yBkgn)
        return self.yBkgn

    def linear(self,slope,x,b):
        return (slope*x)+b


    def better_shirley(self,x,y, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, tot_area):
        #Not sure where this code came from. May need to cite it later. It has similar structure to the shirley formula described in Herrera's paper
        self.y = y
        self.x = x

        E = x
        J = y


        def integralOne(E, J, B, E1=0, E2=-1):
            integral = []
            if E2 < 0:
                E2 = len(J) + E2
            integral = sum([J[n] - B[n] for n in range(E,E2)])
            return integral

        def integralTwo(E, I, B, E1=0, E2=-1):
            integral = []
            if E2 < 0:
                E2 = len(I) + E2
            integral = sum([I[n] - B[n] for n in range(E1,E2)])
            return integral

        def getBn(E,I,B,E1=0,E2=-1):
            I2 = I[E2]
            I1 = I[E1]
            value = I2 + (I1 - I2)/(integralTwo(E,I,B,E1,E2))*integralOne(E,I,B,E1,E2)
            return value

        def iterateOnce(I,B,E1=0,E2=-1):
            b = [getBn(E,I,B,E1,E2) for E in range(len(I))]
            return b

        Bn = [0 for i in range(len(J))]
        Bn = iterateOnce(J,Bn)
        for i in range(6): #how many iterations it's doing
            B_temp = Bn
            Bn = iterateOnce(J,Bn)
            B_diff = [Bn[j] - B_temp[j] for j in range(len(Bn))] #Could make a check to see if the iterations are getting better. Usually little difference after 7 iterations

        self.yBkgn = Bn


        return self.yBkgn








    '''
    def SVSC_shirley(self,x,y, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split):

        self.y = y
        self.x = x

        #Should probably declare these items outside of each background type
        numPeaksUsed = 1
        maxPeaks = 1
        numberBackgrounds = 1
        ma = maxPeaks + numberBackgrounds #will need to change to make it able to use multiple peaks/backgrounds
        numP = len(self.y)
        a =[0]*len(self.x) #dont know if we need a anymore?

        #voigt = peak.voigt
        funcs = y #len = numP -1 #recover initial peak curve --> How to get y points of just one peak??? Use BE range + some delta
        backgroundFromPeakShirleyFix = [0]*(len(self.y)-1) #not sure why its one less than the number of points
        SVSC_bkgn = backgroundFromPeakShirleyFix #easier to write --> original name comes from aanalyzer code

        a_old = 0.3
        a_new = 0.5 #are these initial values too large?
        old_fit = 10000
        best_fit = funcs #setting initial best fit --> just equal to y originally
        SVSC_diff = 1
        while a_new >= 0: #Iterates until a = 0, but keeps track of std of background to voigt fit. Need to find a better way for the GA to optimize a
            i = 1
            for i in range(maxPeaks):#calculates background for each peak then iterates
                #a_ratio is some parameter ratio --> I think it is the ratio of one parameter of different correlated peaks, unsure as to which parameter is being correlated
                a_ratio_b4 = a_old #Right now these are just random --> real code: a[ peakShirleyma[ peakShirelyCorrTo] ] / a[ mama[peakShirelyCorrTo] ]
                a_ratio_after = a_new #defined on line 15233 in PUnit1 --> Values are a[] before and after lfitmod is called
                peakShirleyBackground = 0.8*a_ratio_b4 + 0.2*a_ratio_after #I think this is supposed to be the scattering factor? Now sure how it is optimized
                #Maybe for now we should treat peakShirleyBackground as the scattering factor?

                for j in np.arange(numP -2, 0, -1):
                    SVSC_bkgn[j-1] = self.y[j-1]*-(self.x[j+1]-self.x[j])*peakShirleyBackground + SVSC_bkgn[j] #isnt this just what we already had but now with a wider range?
                    funcs[j] += SVSC_bkgn[j-1]
                #should write array in here to store each peak curves background --> will sum these up later
                i +=1

                iteration_diff = np.subtract(voigt, funcs) #need to change voigt to whatever the curve fit y array is
                new_fit = np.std(iteration_diff)
                new_fit_array = funcs
                if new_fit < old_fit:
                    old_fit = new_fit
                    best_fit = new_fit_array
                a_old = a_new
                a_new -= 0.01 #slow decrease for now --> NEED TO FIND BETTER WAY TO OPTIMIZE a_new
                #lfitmod caluculated here --> Calcualtes parameters between iterations: This is what makes the background active
                #Should we call class Peak here to recalculate the fit with the new background? Active curve fitting

        funcs = best_fit
        #return funcs #Not sure how we are calling this (self.yBKgn?)
        self.yBkgn = funcs

        return self.yBkgn
    '''





    def shirley_Sherwood(self,x,y, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, tot_area):
        #too lazy to find all the x and ys and get rid of the self
        self.y = y
        self.x = x

        #need to make these parameters outside background functions
        numPeaksUsed = 1
        maxPeaks = 1
        numberBackgrounds = 1
        ma = maxPeaks + numberBackgrounds
        useIntegralBkgn=True
        numP = [0]*len(self.y) #we are using this as an array right now but it should just be the number of data points
        a =[0]*len(self.x)


        def iterations(self,x,y):
            numPeaksUsed = 1
            maxPeaks = 1
            numberBackgrounds = 1
            ma = maxPeaks + numberBackgrounds
            useIntegralBkgn=True
            numP = [0]*len(self.y) #we are using this as an array right now but it should just be the number of data points
            a =[0]*len(self.x)
            #need this to find the correct data points in which the bakcground will be removed
            numPointsAroundBackgroundLimitsLocal = 5
            nRightLocal = numPointsAroundBackgroundLimitsLocal // 2
            nLeftLocal = numPointsAroundBackgroundLimitsLocal // 2

            yRightLocal = 0
            yLeftLocal = 0

            for j in range(-(numPointsAroundBackgroundLimitsLocal // 2), numPointsAroundBackgroundLimitsLocal // 2 + 1):
                #yRightLocal += datos[dataNumber].ModifiedCurve.y[nRightLocal + j]
                #yLeftLocal += datos[dataNumber].ModifiedCurve.y[nLeftLocal + j]
                yRightLocal += self.y[len(self.y) - nRightLocal-1 + j]
                yLeftLocal += self.y[nLeftLocal + j]

            yLeftLocal /= (numPointsAroundBackgroundLimitsLocal // 2) * 2 + 1
            yRightLocal /= (numPointsAroundBackgroundLimitsLocal // 2) * 2 + 1



            nLeft = nLeftLocal
            nRight = len(x)-nRightLocal
            global yRight
            yRight = yRightLocal
            yLeft = yLeftLocal


            iterationsIntegralBkgn = 6
            BkgdCurve = []
            funcs = numP
            #funcs = np.zeros((ma, numP))

            if useIntegralBkgn: #from Aanalyzer code line 14867
                ma += 1 #Need to add in array to store each peak background peak[ma] --> sum in other backgrounds? peak[ma] += funcs...
                #funcs[ma][numP] = 0
                #print("K is " + str(self.k))
                for j in range(nRight-1, -1, -1):
                    #funcs[ma][j] = (self.y[j] - yRight[j]) * (self.x[j+1] - self.x[j]) + funcs[ma][j+1]
                    #print(self.y[j] - yRight)
                    funcs[j] = (self.y[j] - yRight) *self.k* -(self.x[j+1] - self.x[j]) + funcs[j+1] #assumes x is in KE, not sure if that changes anything

                '''
                for j in range(0, nLeft):
                    #funcs[ma][j] = funcs[ma][nLeft]
                    funcs[j] = yLeft-yRight
                '''
                integralma = ma

            '''
            #iterates shirley background
            if useIntegralBkgn: #from Aanalyzer code line 15140
                for l in range(iterationsIntegralBkgn):
                    for j in range(nRight-1, nLeft, -1):
                        #funcs[integralma][j] = (self.y[j] - yRight[j] - a[integralma] * funcs[integralma][j]) * (self.x[j+1] - self.x[j]) + funcs[integralma][j+1]
                        funcs[j] = (self.y[j] - yRight - funcs[j]) * (self.x[j+1] - self.x[j]) + funcs[j+1]
                    for j in range(1, nLeft):
                        #funcs[integralma][j] = funcs[integralma][nLeft]
                        funcs[j] = funcs[nLeft]
                        #calls lfitmod here -->calculates chisq and deletes all parameters

                    l += 1
            '''
            return funcs
        for i in range(1): #How many iterations it is performing
            funcs = iterations(self,x,y)

        self.yBkgn = funcs

        ''' Old built in baseline (bad)
        for i in range(len(self.yBkgn)):
            self.yBkgn[i] += yRight

        return self.yBkgn
        '''
    '''
    Just barely started on peak shirley, commented out so it wont cause a compilation error
    def peak_shirley(self,x,y,peak):
        peak.getY
    '''
    @nb.jit
    def Tougaard3Param(self,x,y, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, tot_area):
        funcs = np.zeros(len(y))
        numP = len(y) #we are using this as an array right now but it should just be the number of data points

        #need this to find the correct data points in which the bakcground will be removed
        numPointsAroundBackgroundLimitsLocal = 5
        nRightLocal = numPointsAroundBackgroundLimitsLocal // 2
        nLeftLocal = numPointsAroundBackgroundLimitsLocal // 2

        yRightLocal = 0
        yLeftLocal = 0

        for j in range(-(numPointsAroundBackgroundLimitsLocal // 2), numPointsAroundBackgroundLimitsLocal // 2 + 1):
            yRightLocal += y[len(y) - nRightLocal-1 + j]
            yLeftLocal += y[nLeftLocal + j]

        yLeftLocal /= (numPointsAroundBackgroundLimitsLocal // 2) * 2 + 1
        yRightLocal /= (numPointsAroundBackgroundLimitsLocal // 2) * 2 + 1

        nRightLocal = len(x)-nRightLocal

        CTou3 = 1000
        DTou3 = 13300

        for i in np.arange(nLeftLocal, nRightLocal):
            integralGar = 0
            for j in np.arange(i, nRightLocal):
                energyDif = -(x[j] - x[i])
                integralGar += (y[j] - yRightLocal)*energyDif/(pow(CTou3-pow(energyDif,2),2) + (DTou3*pow(energyDif,2))*-(x[j+1]-x[j]))
            funcs[i] = integralGar

        BTou3 = (yLeftLocal - yRightLocal)/funcs[nLeftLocal]
        for i in np.arange(nLeftLocal, nRightLocal):
            funcs[i] *= BTou3
            funcs[i] += yRightLocal

        for i in np.arange(0, nLeftLocal):
            funcs[i] = yLeftLocal

        for i in np.arange(nRightLocal, numP):
            funcs[i] = yRightLocal
            self.yBkgn = funcs

        return self.yBkgn





    @nb.jit
    def Tougaard2Param(self,x,y, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, tot_area): #Added in but does not work right now
        funcs = np.zeros(len(y))
        numP = len(y) #we are using this as an array right now but it should just be the number of data points

        #need this to find the correct data points in which the bakcground will be removed
        numPointsAroundBackgroundLimitsLocal = 5
        nRightLocal = numPointsAroundBackgroundLimitsLocal // 2
        nLeftLocal = numPointsAroundBackgroundLimitsLocal // 2

        yRightLocal = 0
        yLeftLocal = 0

        for j in range(-(numPointsAroundBackgroundLimitsLocal // 2), numPointsAroundBackgroundLimitsLocal // 2 + 1):
            yRightLocal += y[len(y) - nRightLocal-1 + j]
            yLeftLocal += y[nLeftLocal + j]

        yLeftLocal /= (numPointsAroundBackgroundLimitsLocal // 2) * 2 + 1
        yRightLocal /= (numPointsAroundBackgroundLimitsLocal // 2) * 2 + 1

        nRightLocal = len(x)-nRightLocal

        CTou2 = 1643 #Taken from Aanalyzer. Given value in the literature.



        for i in np.arange(nLeftLocal, nRightLocal):
            integralGar = 0
            for j in np.arange(i, nRightLocal):
                energyDif = -(x[j] - x[i])
                integralGar += (y[j] - yRightLocal)*energyDif/CTou2-pow(energyDif,2)/CTou2-pow(energyDif,2)*(x[j+1] - x[j])
            funcs[i] = integralGar

        BTou2 = (yLeftLocal - yRightLocal)/funcs[nLeftLocal] #Should be close to 3000 for most metals
        for i in np.arange(nLeftLocal, nRightLocal):
            funcs[i] *= BTou2
            funcs[i] += yRightLocal

        for i in np.arange(0, nLeftLocal):
            funcs[i] = yLeftLocal

        for i in np.arange(nRightLocal, numP):
            funcs[i] = yRightLocal
            self.yBkgn = funcs

        return self.yBkgn






    def shirley_bkgn_again(self,x,y, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, tot_area):

       
        #Dont want to be dependent upon I1 and I2:
        #Change to be more like:
        #
        #
        #S(E) = baseline_shirley + k1*(intensity integral peak 1)*Area peak 1 + k2*(intensity integral peak 2)*Area peak 2 .... for each peak
        #
        #Each peak is given a unique scattering factor valued between 0 and 1
        self.y = y
        self.x = x
  
        #integral = [0]*len(self.y) #we are using this as an array right now but it should just be the number of data points\
        #integral = [self.backgroundShirley]*len(y)

       
     
        #Phase 4.2: container-mode throwaway (no ranges, NO RNG DRAWS) - the
        #drawn values were never used; the shape function below receives all
        #parameters explicitly. This removes RNG consumption from fitness
        #evaluation, making parallel evaluation trajectory-identical to
        #serial (goldens rebaselined; see CHANGELOG).
        peak_class = peak('', self.peakType)
        #baseline_shirley = [y[-1]]*len(y) #First value is just last data point in y array aka the rightside of data
        #self.bindingEnergy = peak_class.bindingEnergy
      
        
        if(self.peakType.lower() == "voigt"):


            self.yValues = peak_class.voigtFunc(x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split) #These are not the same as the ones we get from peak.peakFunc --> Object attribute error

        elif(self.peakType.lower() == "gaussian"):
            self.yValues = peak_class.gaussFunc(x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split)
        elif(self.peakType.lower() == "lorentzian"):
            self.yValues = peak_class.lorentzFunc(x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split)
        elif(self.peakType.lower() == "double lorentzian"):
            
            self.yValues = peak_class.doubleLorentzFunc(x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split)

        elif(self.peakType.lower() == "doniach-sunjic"):
            self.yValues = peak_class.doniachSunjicFunc(x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split)

        else:
            print("Error assigning peak type")
            print("Peaktype found is: " + str(self.peakType))
            sys.exit(1)
        
        baseline_shirley = [min(self.yValues)]*len(y)
        #baseline_shirley = [min(self.y)]*len(y) #Using the smallest y value as the basis for the right side of the data
        #yValues = tot_area
        yValues = self.yValues
        

        '''
        for i in range(len(x)):
            print(x[i], self.yValues[i])
        '''
    
        def areaA(baseline_shirley):
            numPointsAroundBackgroundLimitsLocal = 5
            nRightLocal = numPointsAroundBackgroundLimitsLocal // 2


            yRightLocal = 0


            for j in range(-(numPointsAroundBackgroundLimitsLocal // 2), numPointsAroundBackgroundLimitsLocal // 2 + 1):
                yRightLocal += self.y[len(self.y) - nRightLocal-1 + j]
            #Upper bound of integral
            yRightLocal /= (numPointsAroundBackgroundLimitsLocal // 2) * 2 + 1
            nRight = len(x)-nRightLocal
            #Phase 4.2 (verified equivalent): vectorized trapezoid segments
            #with np.cumsum, which accumulates sequentially in the same
            #descending-i order as the former scalar loop.
            areaA = np.zeros(len(y))
            area_A = np.zeros(len(y))
            _yv = np.asarray(self.yValues, dtype=np.float64)
            _xv = np.asarray(x, dtype=np.float64)
            _seg = ((_yv[:nRight] + _yv[1:nRight+1]) / 2) * (_xv[1:nRight+1] - _xv[:nRight])
            areaA[:nRight] = _seg
            _acc = np.cumsum(_seg[::-1])
            area_A[:nRight] = np.concatenate(([0.0], _acc[:-1]))[::-1]
            return abs(area_A)

        total_area = areaA(baseline_shirley)

        #baseline_shirley = [min(self.y)]*len(y)

        def iterations(self,x,baseline_shirley, k_val, yValues):
            
            #need this to find the correct data points in which the bakcground will be removed
            baseline_shirley = [min(yValues)]*len(y)
            numPointsAroundBackgroundLimitsLocal = 5
            nRightLocal = numPointsAroundBackgroundLimitsLocal // 2


            yRightLocal = 0


            for j in range(-(numPointsAroundBackgroundLimitsLocal // 2), numPointsAroundBackgroundLimitsLocal // 2 + 1):
                yRightLocal += self.y[len(self.y) - nRightLocal-1 + j]
            #Upper bound of integral
            yRightLocal /= (numPointsAroundBackgroundLimitsLocal // 2) * 2 + 1
            nRight = len(x)-nRightLocal
            yRight = yRightLocal



            #N = 30
            #y_right = y[-N:]
            #y_right_avg = sum(y_right)/N
            #base = y_right_avg - tot_area[-1]
            #integral = [base]*len(y)
            integral = baseline_shirley
            new_integral = baseline_shirley
            
            #new_integral = [base]*len(y)
            deltaX = abs(x[1] - x[2]) #Data stepsize --> this is assumed that the x data is stepwise and equal throughout the whole data range
            #Integral Calculation

            #Is nRight the correct value?           
            #Phase 4.2 (verified equivalent): elementwise update vectorized;
            #each slot was written once from its own original value.
            new_integral = np.asarray(new_integral, dtype=np.float64)
            new_integral[:nRight] = np.asarray(integral[:nRight], dtype=np.float64) + k_val*np.asarray(total_area[:nRight], dtype=np.float64)
                
            #print(np.array(new_integral).sum(), np.array(integral).sum())
            #return integral
            return new_integral

        self.yBkgn = iterations(self,x, baseline_shirley, self.k, yValues)

        '''
        for i in range(1):
       
            self.yBkgn = iterations(self,x, baseline_shirley, self.k, yValues)
            
            area_before = total_area
            total_area = areaA(self.yBkgn)
           
            #bkgn_diff = self.yBkgn[-1] - self.y[-1]
            #baseline_shirley = [min(self.y) - bkgn_diff]*len(y)
        

            
            chance = 20
            if random.random()*100 < chance:
                # for j in range left-right
                #
                # k = y_left - y_right / integral[j] + k_val*total_area[j] - new_integral[j]
                #
                diff = 0
                numPointsAroundBackgroundLimitsLocal = 5
                for j in range(-(numPointsAroundBackgroundLimitsLocal // 2), numPointsAroundBackgroundLimitsLocal // 2 + 1):
                    diff += area_before[j] - self.yBkgn[j] #self.k*total_area[j]
                    #diff += total_area[j] - self.yBkgn[j]

                self.k = (y[0]-y[-1])/(diff)
                if self.k < 0: #No negative k
                    self.k = 0.0
         
            self.k = (y[0]-y[-1])/((area_before[0] - self.k*total_area[0]))
            if self.k < 0: #No negative k
                self.k = 0.0
            
            #self.k = (y[0]-y[-1])/((total_area[0] - self.k*total_area[0]))
        '''


        return self.yBkgn
    




    #Regular Shirley: Does not have different k values for each peak.


    def shirley(self,x,y, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, tot_area):
        
        #Dont want to be dependent upon I1 and I2:
        #Change to be more like:
        #
        #
        #S(E) = baseline_shirley + k1*(intensity integral peak 1)*Area peak 1 + k2*(intensity integral peak 2)*Area peak 2 .... for each peak
        #
        #Each peak is given a unique scattering factor valued between 0 and 1
        self.y = y
        self.x = x
  
        #integral = [0]*len(self.y) #we are using this as an array right now but it should just be the number of data points\
        #integral = [self.backgroundShirley]*len(y)

        
        #Phase 4.2: container-mode throwaway (no ranges, NO RNG DRAWS) - the
        #drawn values were never used; the shape function below receives all
        #parameters explicitly. This removes RNG consumption from fitness
        #evaluation, making parallel evaluation trajectory-identical to
        #serial (goldens rebaselined; see CHANGELOG).
        peak_class = peak('', self.peakType)
        #baseline_shirley = [y[-1]]*len(y) #First value is just last data point in y array aka the rightside of data
        #self.bindingEnergy = peak_class.bindingEnergy
      

        if(self.peakType.lower() == "voigt"):


            self.yValues = peak_class.voigtFunc(x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split) #These are not the same as the ones we get from peak.peakFunc --> Object attribute error

        elif(self.peakType.lower() == "gaussian"):
            self.yValues = peak_class.gaussFunc(x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split)
        elif(self.peakType.lower() == "lorentzian"):
            self.yValues = peak_class.lorentzFunc(x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split)
        elif(self.peakType.lower() == "double lorentzian"):
            
            self.yValues = peak_class.doubleLorentzFunc(x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split)

        elif(self.peakType.lower() == "doniach-sunjic"):
            self.yValues = peak_class.doniachSunjicFunc(x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split)

        else:
            print("Error assigning peak type")
            print("Peaktype found is: " + str(self.peakType))
            sys.exit(1)
        baseline_shirley = [min(self.yValues)]*len(y)
        #baseline_shirley = [min(self.y)]*len(y) #Using the smallest y value as the basis for the right side of the data
        #baseline_shirley = [0]*len(y)
        #yValues = self.yValues
        

        yValues = tot_area
        self.yValues = tot_area

        '''
        for i in range(len(x)):
            print(x[i], self.yValues[i])
        '''


        def areaA(baseline_shirley):

            numPointsAroundBackgroundLimitsLocal = 5
            nRightLocal = numPointsAroundBackgroundLimitsLocal // 2


            yRightLocal = 0


            for j in range(-(numPointsAroundBackgroundLimitsLocal // 2), numPointsAroundBackgroundLimitsLocal // 2 + 1):
                yRightLocal += self.y[len(self.y) - nRightLocal-1 + j]
            #Upper bound of integral
            yRightLocal /= (numPointsAroundBackgroundLimitsLocal // 2) * 2 + 1
            nRight = len(x)-nRightLocal
            #Phase 4.2 (verified equivalent): vectorized trapezoid segments
            #with np.cumsum, which accumulates sequentially in the same
            #descending-i order as the former scalar loop.
            areaA = np.zeros(len(y))
            area_A = np.zeros(len(y))
            _yv = np.asarray(self.yValues, dtype=np.float64)
            _xv = np.asarray(x, dtype=np.float64)
            _seg = ((_yv[:nRight] + _yv[1:nRight+1]) / 2) * (_xv[1:nRight+1] - _xv[:nRight])
            areaA[:nRight] = _seg
            _acc = np.cumsum(_seg[::-1])
            area_A[:nRight] = np.concatenate(([0.0], _acc[:-1]))[::-1]
            return abs(area_A)

        total_area = areaA(baseline_shirley)
        


     
        def iterations(self,x,baseline_shirley, k_val, yValues):

            #need this to find the correct data points in which the bakcground will be removed
            baseline_shirley = [min(yValues)]*len(y)
            numPointsAroundBackgroundLimitsLocal = 5
            nRightLocal = numPointsAroundBackgroundLimitsLocal // 2


            yRightLocal = 0

          
            for j in range(-(numPointsAroundBackgroundLimitsLocal // 2), numPointsAroundBackgroundLimitsLocal // 2 + 1):
                yRightLocal += self.y[len(self.y) - nRightLocal-1 + j]
            #Upper bound of integral
            yRightLocal /= (numPointsAroundBackgroundLimitsLocal // 2) * 2 + 1
            nRight = len(x)-nRightLocal
            yRight = yRightLocal


            N = 30
            y_right = y[-N:]
            y_right_avg = sum(y_right)/N
            #integral = baseline_shirley
            base = y_right_avg - tot_area[-1] #This garentees that the background passes through the points on the right (low BE) side of data. Should this be allowed to change/be picked by the algorithm to allow for lower values? --> Really low baseline? 
            integral = [base]*len(y)
            

            new_integral = [base]*len(y)

            #new_integral = baseline_shirley
            deltaX = abs(x[1] - x[2]) #Data stepsize --> this is assumed that the x data is stepwise and equal throughout the whole data range
            #Integral Calculation

            #Phase 4.2 (verified equivalent): elementwise update vectorized;
            #each slot was written once from its own original value.
            new_integral = np.asarray(new_integral, dtype=np.float64)
            new_integral[:nRight] = np.asarray(integral[:nRight], dtype=np.float64) + k_val*np.asarray(total_area[:nRight], dtype=np.float64)
                #print(k_val*total_area[j])
                #print(k_val*total_area[j])
            #print(np.array(new_integral).sum(), np.array(integral).sum())
            #return integral
          
            return new_integral

        

        for i in range(1):
       

            
            self.yBkgn = iterations(self,x, baseline_shirley, self.k, yValues)
          
            '''
            area_before = total_area
            total_area = areaA(self.yBkgn)
            
            chance = 20
            if random.random()*100 < chance:
                # for j in range left-right
                #
                # k = y_left - y_right / integral[j] + k_val*total_area[j] - new_integral[j]
                #
                diff = 0
                numPointsAroundBackgroundLimitsLocal = 5
                for j in range(-(numPointsAroundBackgroundLimitsLocal // 2), numPointsAroundBackgroundLimitsLocal // 2 + 1):
                    diff += area_before[j] - self.yBkgn[j] #self.k*total_area[j]
                    #diff += total_area[j] - self.yBkgn[j]

                self.k = (y[0]-y[-1])/(diff)
                if self.k < 0: #No negative k
                    self.k = 0.0
                
               
                #self.k = (y[0]-y[-1])/((total_area[0] - self.k*total_area[0]))
            '''

        return self.yBkgn

    #----------------- Phase 3d: analysis value extractors -----------------
    # Ported from the gui fork (get_2Param_vals / get_3Param_vals) with the
    # peakType plumbing dropped; used by Individual.get_analysis_params to
    # report Tougaard B values. Constants are the Aanalyzer defaults.

    def get_2Param_vals(self, x, y, tot_peak_fit):
        funcs = np.zeros(len(y))
        numPointsAroundBackgroundLimitsLocal = 5
        nRightLocal = numPointsAroundBackgroundLimitsLocal // 2
        nLeftLocal = numPointsAroundBackgroundLimitsLocal // 2

        yRightLocal = 0
        yLeftLocal = 0

        for j in range(-(numPointsAroundBackgroundLimitsLocal // 2), numPointsAroundBackgroundLimitsLocal // 2 + 1):
            yRightLocal += y[len(y) - nRightLocal-1 + j]
            yLeftLocal += y[nLeftLocal + j]

        yLeftLocal /= (numPointsAroundBackgroundLimitsLocal // 2) * 2 + 1
        yRightLocal /= (numPointsAroundBackgroundLimitsLocal // 2) * 2 + 1

        nRightLocal = len(x)-nRightLocal

        CTou2 = 1643 #Taken from Aanalyzer

        for i in np.arange(nLeftLocal, nRightLocal):
            integralGar = 0
            for j in np.arange(i, nRightLocal):
                energyDif = -(x[j] - x[i])
                integralGar += (y[j] - yRightLocal)*energyDif/CTou2-pow(energyDif,2)/CTou2-pow(energyDif,2)*(x[j+1] - x[j])
            funcs[i] = integralGar

        BTou2 = (yLeftLocal - yRightLocal)/funcs[nLeftLocal]
        self.BTou2 = BTou2
        return BTou2

    def get_3Param_vals(self, x, y, tot_peak_fit):
        funcs = np.zeros(len(y))
        numPointsAroundBackgroundLimitsLocal = 5
        nRightLocal = numPointsAroundBackgroundLimitsLocal // 2
        nLeftLocal = numPointsAroundBackgroundLimitsLocal // 2

        yRightLocal = 0
        yLeftLocal = 0

        CTou3 = 1000
        DTou3 = 13300

        for j in range(-(numPointsAroundBackgroundLimitsLocal // 2), numPointsAroundBackgroundLimitsLocal // 2 + 1):
            yRightLocal += y[len(y) - nRightLocal-1 + j]
            yLeftLocal += y[nLeftLocal + j]

        yLeftLocal /= (numPointsAroundBackgroundLimitsLocal // 2) * 2 + 1
        yRightLocal /= (numPointsAroundBackgroundLimitsLocal // 2) * 2 + 1

        nRightLocal = len(x)-nRightLocal

        for i in np.arange(nLeftLocal, nRightLocal):
            integralGar = 0
            for j in np.arange(i, nRightLocal):
                energyDif = -(x[j] - x[i])
                integralGar += (y[j] - yRightLocal)*energyDif/(pow(CTou3-pow(energyDif,2),2) + (DTou3*pow(energyDif,2))*-(x[j+1]-x[j]))
            funcs[i] = integralGar

        BTou3 = (yLeftLocal - yRightLocal)/funcs[nLeftLocal]
        self.BTou3 = BTou3
        return BTou3
