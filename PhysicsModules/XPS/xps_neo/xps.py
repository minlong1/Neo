import enum
from PhysicsModules.XPS.xps_neo.helper import *
from PhysicsModules.XPS.xps_neo.import_lib import *
from PhysicsModules.XPS.xps_neo.input_arg import parse_args, load_file_dict
from PhysicsModules.XPS.xps_neo.ini_parser import load_config
from PhysicsModules.XPS.xps_neo.loss import compute_loss
#from pathObj import OliverPharr
#from .individual import Individual
#from pathrange import Pathrange_limits  #may be deletable/ built for XAFS
#from xps_neo_data import XPS_Data #Do we need this still?

import numpy as np
from PhysicsModules.XPS.xps_neo.xps_individual import Individual
from PhysicsModules.XPS.xps_neo.xps_fit import peak,background
from PhysicsModules.XPS.xps_neo.xps_data import xps_data

from copy import deepcopy #fixes bug at line 70ish with deepcopy
"""
Author: Alaina Humiston, Miu Lun Lau
"""


class XPS_GA:
    """XPS Neo's self-contained GA/DE driver: owns the heterogeneous-genome
    population (see the module's own `individual.get_params()` layout —
    floats interleaved with peak/background type-name strings), the
    add/remove-peak logic, and the generation loop (`run`). Its methods
    deliberately read fit configuration as bare globals of this module
    (set once via `globals().update(config)` in `main()`), a documented
    transitional seam carried over from the source project rather than an
    oversight — see PhysicsModules/XPS/README.md "Why this module doesn't
    use Solvers" for why it doesn't route through Solvers.core."""

    def initialize_params(self,verbose = False):
        """
        Initialize Parameters
        """
        # print("Initialize Parameters")
        print("Initializing Params")
        self.intervalK = 0.05
        self.numGenSinceImproved = 0
        self.tol = np.finfo(np.float64).resolution
        

    def initialize_variable(self):
        """
        Initialize variables
        """
        print("Initializing Variables")
        self.genNum = 0
        self.nChild = 4
        self.globBestFit = [0,0]
        self.currBestFit = [0,0]
        self.bestDiff = 9999e11
        self.bestBest = 999999999e11
        self.diffCounter = 0
        self.element = element
        self.photoelectronLine = photoelectronLine
        self.pathDictionary = {}
        self.data_file = data_file
        #self.data_cutoff = data_cutoff
        # Paths
        self.npaths = npaths
        #self.fits = fits

        # Populations
        self.npops = size_population
        self.ngen = number_of_generation
        self.steady_state = steady_state

        #Alter 
        self.gen_alt = gen_alt

        # Mutation Parameters
        self.mut_opt = mutated_options
        self.mut_chance = chance_of_mutation
        # self.mut_chance_e0 = chance_of_mutation_e0

        # Crosover Parameters
        self.n_bestsam = int(best_sample*self.npops*(0.01))
        self.n_lucksam = int(lucky_few*self.npops*(0.01))

        # Time related
        self.time = False
        self.tt = 0

        # DE related:
        self.F = 0.5
        self.cR = 0.3

    def initialize_file_path(self,i=0):
        """
        Initialize file paths for each of the file first
        """
        print("initializing file paths")
        self.base = os.getcwd()
        self.output_path = os.path.join(self.base,output_file)
        self.check_output_file(self.output_path)
        self.log_path = os.path.splitext(copy.deepcopy(self.output_path))[0] + ".log"
        self.check_if_exists(self.log_path)

        # Initialize logger (named, not the root logger, so xps_neo does
        # not hijack the logging of embedding applications - Phase 5)
        self.logger = logging.getLogger('xps_neo')
        self.logger.propagate = False
        # Delete handler
        self.logger.handlers=[]
        file_handler = logging.FileHandler(self.log_path,mode='a+',encoding='utf-8')
        stdout_handler = logging.StreamHandler(sys.stdout)

        formatter = logging.Formatter('%(message)s')
        file_handler.setFormatter(formatter)
        stdout_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)
        self.logger.addHandler(stdout_handler)

        self.logger.setLevel(logging.DEBUG if getattr(self, "verbose", False) else logging.INFO)
        self.logger.info(banner())

    def check_if_exists(self,path_file):
        """
        Check if the directory exists
        """
        if os.path.exists(path_file):
            os.remove(path_file)
        # Make Directory when its missing
        path = pathlib.Path(path_file)
        path.parent.mkdir(parents=True, exist_ok=True)

    def check_output_file(self,file):
        """
        check if the output file for each of the file
        """
        file_base= os.path.splitext(file)[0]
        self.check_if_exists(file)
        self.file = file

        self.file_initial = open(self.output_path,"a+")
        self.file_initial.write("Gen,TPS,FITTNESS,CURRFIT,CURRIND,BESTFIT,BESTIND\n")  # writing header
        self.file_initial.close()

        file_data = os.path.splitext(file)[0] + '_data.csv'
        self.check_if_exists(file_data)
        self.file_data = file_data


    def initialize_range(self,scale_var, x_offset, y_offset, data_KE, data_XES,i=0,BestIndi=None):
       
        print("initializing ranges")
        """
        Initalize range

        To Do list:
            Initalize range will be difference for each paths depend if the run are
            in series, therefore the ranges will self-adjust
        """
        
  
        photoline_select = self.photoelectronLine
        scale_var = scale_data
        x_offset = x_offset
        y_offset = y_offset
        self.data_obj = xps_data(self.data_file,skipLn = skipLn, x_offset=x_offset, y_offset=y_offset)
        #HOW TO CALL SCALE_VAR AND GET CORRECT VALUE IN THIS FUNCTION OR INIT?
        self.x_slice = self.data_obj.get_x()
        self.y_slice = self.data_obj.get_y(scale_var)
        self.x_array = self.x_slice
        self.y_array = self.y_slice
       
        yAvg = 0
        yTot = 0
        j=0
        for i,yVal in enumerate(self.data_obj.get_y(scale_var)[-10:]): #[-10:] gets last 10 items in array.
            yTot += yVal
            j=i
        yAvg = yTot/(j+1)

        #baseline_range[1] = yAvg
        #[0] = min range, [1] = max range, [2] = delta range
        #baseline_range = [0,yAvg,1]
        #amp_range[1] = max(self.data_obj.get_y()) #From data tesing it seems that having this in actually helps a lot with finding parameters
        

        #Added in a lot of ranges and if statements to help limit the background to the data. Could cause some errors. 
        N = 30
        y_left = self.y_array[:N]
        y_right = self.y_array[-N:]
        x_left = self.x_array[:N]
        x_right = self.x_array[-N:]

        self.y_left_avg = sum(y_left)/N
        x_left_avg = sum(x_left)/N
        self.y_right_avg = sum(y_right)/N
        x_right_avg = sum(x_right)/N
        if self.y_left_avg > self.y_right_avg:
            #Should we scale baseline so that the last point (lowest BE/highest KE) is subtratced off the y values (goes to zero) and set highest possible baseline value to the left side of the data?
            #baseline_range = [(self.y_right_avg-np.round(self.y_right_avg*0.05, 2)),self.y_right_avg+np.round(self.y_right_avg*0.05, 2),0.1] #Changed minimum to allow for background lower. Now allows for 20% below avg. and 5 % above
            baseline_range = [(self.y_right_avg-np.round(self.y_right_avg*0.25, 2)),self.y_right_avg,0.1]
            #baseline_range = [-1000.0,10.0,1]
        else:
            #baseline_range = [(self.y_left_avg-np.round(self.y_left_avg*0.05, 2)),self.y_left_avg+np.round(self.y_left_avg*0.05, 2),0.1] #KE option
            baseline_range = [(self.y_left_avg-np.round(self.y_left_avg*0.25, 2)),self.y_left_avg,0.1]
            #baseline_range = [1000,self.y_right_avg+np.round(self.y_right_avg*0.01, 2),1]

        if 'Linear' in background_type: #set baseline to zero
            baseline_range = [-0.001, 0.001, 0.001]
        elif '3-Param Tougaard'in background_type: #set baseline to zero
            baseline_range = [-0.001, 0.001, 0.001]
        elif '2-Param Tougaard'in background_type: #set baseline to zero
            baseline_range = [-0.001, 0.001, 0.001]


        '''
        if 'SVSC' and 'Shirley'in background_type:
            background_type.remove('Shirley') #Only need one Shirley background
        '''


        if 'SVSC' in background_type: #Code throws an error if Shirley isn't last in background array. Issue is found in gui/xps_individual line 134 
            background_type.remove('SVSC')
            background_type.append('SVSC')
        if 'Shirley' in background_type: #Code throws an error if Shirley isn't last in background array. Issue is found in gui/xps_individual line 134 
            background_type.remove('Shirley')
            if 'Baseline' in background_type:
                background_type.remove('Baseline')
            background_type.append('Shirley')


        y_slope = (self.y_right_avg - self.y_left_avg)/abs(x_right_avg - x_left_avg)
       


        if y_slope < 0:
           
            slope_range[0] = y_slope*1.1
            slope_range[1] = y_slope*0.9
            if 'SVSC' or 'Shirley' in background_type:
                slope_range[0] = y_slope*1.0
                slope_range[1] = 0 #y_slope*0.5
                k_range = [0.01, 0.1, 0.001]
                
                background_range[1] = self.y_left_avg
                background_range[0] = self.y_right_avg*0.5 #Smaller if shirley to allow for shirley step function
            else:
                
                k_range = [0, 0.1, 0.001]
                background_range[1] = self.y_left_avg
                background_range[0] = self.y_right_avg

        else:
            slope_range[0] = y_slope*0.9
            slope_range[1] = y_slope*1.1 #Set to the slope of the average of both sides of the data plus a little more for wiggle room
            if 'SVSC' or 'Shirley' in background_type:
                slope_range[1] = y_slope*1.0
                slope_range[0] = 0 #y_slope*0.5
                k_range = [0.01, 0.1, 0.001]
                background_range[0] = self.y_left_avg*0.5
                background_range[1] = self.y_right_avg
            else:
                k_range = [0, 0.1, 0.001]
                background_range[0] = self.y_left_avg
                background_range[1] = self.y_right_avg
       
        '''
        if x_left_avg > x_right_avg:
            background_range[1] = x_left_avg
            background_range[0] = x_right_avg
        else:
            background_range[0] = x_left_avg
            background_range[1] = x_right_avg
        '''
        #Check to make sure that the Shirley background value does not go to zero for some reason. How to better correlate slpoe to k?
        
        '''
        if y_slope >= 10:
            k_range = [0.05, 0.1, 0.001] 
        else: 
            k_range = [0, 0.1, 0.001]
        '''
        


        '''
        if y_slope < 0:
            slope_range[0] = y_slope - 5

        else:
            slope_range[1] = y_slope + 5 #Set to the slope of the average of both sides of the data plus a little more for wiggle room

        background_range[1] = self.y_left_avg + 200
        '''
        
        diff = self.y_left_avg - self.y_right_avg
        #background_shir_range[1] = self.y_right_avg + (self.y_right_avg/100)*25


        #Can change ranges to be within a certain percent range of what the user inputs. Right now it is set at 35% but can increase if needed

        '''
        percent_guess = 35/100

        sigma_range[0] =  -sigma_guess[0]*percent_guess
        sigma_range[1] =  sigma_guess[0]*percent_guess

        fwhm_range[0] = -gamma_guess[0]*percent_guess
        fwhm_range[1] =  gamma_guess[0]*percent_guess

        amp_range[0] =  -self.y_right_avg - amp_guess[0]*percent_guess #subtracts self.y_right_avg because amp guess excludes background (assuming user wont do that correctly)
        amp_range[1] =  -self.y_right_avg + amp_guess[0]*percent_guess
        '''

        #print(amp_guess[0])

        for n in range(len(amp_guess)):
            #amp_guess[n] = amp_guess[n] -round(self.y_right_avg,2) #It is assumed that the user is inputting the amp + baseline so the average of the right data is subtratced to account for this
            if amp_guess[n] < 0:
                amp_guess[n] = abs(amp_guess[n])

        decay_rate = self.y_array[:1]*0.367879441

        for i in np.arange(0, len(self.y_array)): #Finding x value at the decay_rate intensity
            if self.y_array[i] >= decay_rate:
                decay_rate_x = self.x_array[i]
                pass
            else:
                decay_rate_x = self.x_array[i]

        amp_range = []
        asym_range = []
        sos_range = []
        br_range = []
        BE_range = []
        sigma_range = []
        gamma_range = []
        width_guess = []

        """
        for j in range(len(BE_guess)):
            for i in range(len(self.x_array)):
                if self.x_array[i] == round(BE_guess[j], 1):
                    if amp_range_max[j] > self.y_array[i] - amp_guess[j]: #Only adjusts if inputted range is too large
                        amp_range_max[j] = self.y_array[i] - amp_guess[j] #Peaks are only allowed to take on the maximum value of the y value at the given BE
                    print(self.x_array[i], self.y_array[i], BE_guess[j], amp_range_max[j])
        """
        
        
        max_peak = amp_guess[0] #Finding the peak with the largest amplitude and giving it a larger BE range
        #In the future we will need to have the algorithm check the elemental peak assingment here and give it a range based on a database of collected peak information
        peak_location = 0
        for i in amp_guess:
            if i >= max_peak:
                max_peak_location = peak_location
            peak_location += 1
        #BE_range_min[max_peak_location] = -0.5
        #BE_range_max[max_peak_location] = +0.5

        for i in range(len(BE_guess)): #Creating 2D arrays of range values for each peaks BE and amplitude
           
            if amp_limited[i] == True:
                amp_range.append([-0.00001, 0.00001, 0.00001]) #Making the range basically zero if the parameter is 'fixed'
            else:
                amp_range.append([amp_range_min[i], amp_range_max[i], amp_range_delta[i]])



            if BE_limited[i] == True:
                BE_range.append([-0.00001, 0.00001, 0.00001]) #Making the range basically zero if the parameter is 'fixed'
            else:
                BE_range.append([BE_range_min[i], BE_range_max[i], BE_range_delta[i]])



            if sigma_limited[i] == True:
                sigma_range.append([-0.00001, 0.00001, 0.00001]) #Making the range basically zero if the parameter is 'fixed'
            else:
                sigma_range.append([sigma_range_min[i], sigma_range_max[i], sigma_range_delta[i]])



            if gamma_limited[i] == True:
                gamma_range.append([-0.00001, 0.00001, 0.00001]) #Making the range basically zero if the parameter is 'fixed'
            else:
                gamma_range.append([gamma_range_min[i], gamma_range_max[i], gamma_range_delta[i]])

            if asym_limited[i] == True:
                asym_range.append([-0.00001, 0.00001, 0.00001]) #Making the range basically zero if the parameter is 'fixed'
            else:
                asym_range.append([asym_range_min[i], asym_range_max[i], asym_range_delta[i]])

            if sos_limited[i] == True:
                sos_range.append([-0.00001, 0.00001, 0.00001]) #Making the range basically zero if the parameter is 'fixed'
            else:
                sos_range.append([sos_range_min[i], sos_range_max[i], sos_range_delta[i]])
            
            if br_limited[i] == True:
                br_range.append([-0.00001, 0.00001, 0.00001]) #Making the range basically zero if the parameter is 'fixed'
            else:
                br_range.append([br_range_min[i], br_range_max[i], br_range_delta[i]])
            

            
                    
         
    
       
        #For data that comes in KE the spinOrbitSplitting needs to be negative
        #MAY NEED TO EDIT THIS FOR NEW RANGES
        '''
        for i in range(len(is_singlet)):
            if self.x_array[0] < self.x_array[-1]:
                if spinOrbitSplit_guess[i] > 0:
                    spinOrbitSplit_guess[i] = -spinOrbitSplit_guess[i]
        '''

    

        #exp_decay_range[0] = decay_rate_x - decay_rate_x*0.1
        #exp_decay_range[1] = decay_rate_x + decay_rate_x*0.1
        #background_shir_range[0] = y_right_avg - (y_right_avg/100)*25 #I think with the use of multiple peaks this should be allowed to be zero
        #k_range[1] = self.y_left_avg - y_right_avg + 200
        
        self.pars_range = {
            'Binding Energy': BE_range,
            'BE':BE_guess,
            'BE_limited' :BE_limited,
            'BE_correlated' :BE_correlated,
            'BE_correlated_mult' :BE_correlated_mult,
            'Gaussian': sigma_range,
            'Sigma': sigma_guess,
            'sigma_limited' :sigma_limited,
            'sigma_correlated' :sigma_correlated,
            'sigma_correlated_mult' :sigma_correlated_mult,
            'Lorentzian': gamma_range,
            'Gamma': gamma_guess,
            'gamma_limited' :gamma_limited,
            'gamma_correlated' :gamma_correlated,
            'gamma_correlated_mult' :gamma_correlated_mult,
            'Amplitude' : amp_range,
            'Amp': amp_guess,
            'amp_limited' :amp_limited,
            'amp_correlated' :amp_correlated,
            'amp_correlated_mult' :amp_correlated_mult,
            'Asymmetry' : asym_range,
            'Asym': asym_guess,
            'asym_limited' :asym_limited,
            'asym_correlated' : asym_correlated,
            'asym_correlated_mult' :asym_correlated_mult,
            #'Asymmetry' : asymmetry_range,
            'Asymmetry Doniach-Sunjic' : asymmetryDoniach_range,
            'k_range' : k_range,
            'Background' : background_range,
            #'Shirley Background' : background_shir_range,
            'Slope' : slope_range,
            'npeaks' : npaths,
            'baseline' : baseline_range,
            #'branching_ratio_range' : branching_ratio_range,
            #'branching_ratios' : branching_ratio,
            'br_range' : br_range,
            'BR': br_guess,
            'br_limited' :br_limited,
            'br_correlated' : br_correlated,
            'br_correlated_mult' :br_correlated_mult,
            'is_singlet' : is_singlet,
            'is_coster_kronig' : is_coster_kronig,
            'Lorentzian Coster-Kronig' : gamma_CK_range,
            'Gamma Coster-Kronig' : gamma_CK_guess,
            'Asymmetry Range Coster-Kronig' : asym_CK_range,
            'Asymmetry Coster-Kronig' : asym_CK_guess,
            #'spinOrbitSplitting' : spinOrbitSplit_range,
            #'spinOrbitSplit' : spinOrbitSplit_guess, #Are we still using this one??
            'spinOrbitSplitting' : sos_range,
            'SOS': sos_guess,
            'sos_limited' :sos_limited,
            'sos_correlated' : sos_correlated,
            'sos_correlated_mult' :sos_correlated_mult,
            'photoline' : photoline_select
        }
        self.peak_type = peak_type
        self.backgrounds = background_type



    ''' Dont' think it's neccesary, we'll see
    def create_range(self,value,percentage,dt,prec): #-------Where is this called?--------
        """
        Create delta to calculate the ranges
        """
        minus = round(value - percentage*value,prec)
        plus = round(value + percentage*value,prec)
        range = np.arange(minus,plus+dt,dt)
        return range
    '''
   

    def scanResidual(self): #Need to expand function to see if overfitting the data --> Reomve that peak
        """
        Residual-driven peak add/remove decision (part of `peak_add_remove`,
        called from `run` when `data_peak_add` is enabled). See
        PhysicsModules/XPS/README.md "Known limitations" — this path is
        unusable from the CLI, crashing on the first added peak; kept as-is
        per the port's no-fix-while-porting policy.
        """
        #Function to analyze the residual to determine if we need to add/subtract a peak
        
        
        
      
        

        yTotal = np.zeros(len(self.x_array)) 
        self.residual = np.zeros(len(self.x_array))
        residual_noiseless = np.zeros(len(self.x_array))

        yTotal = self.currBestY #Using the best fit Y values of the data 
       
        

        for j in range(len(self.x_array)):
            self.residual[j] = ((self.y_array[j] -yTotal[j])) #Positive means under fitting here and negative means overfitting
       

        #Want to subtract off the avg of the residual from both ends of the data to account for noise in data 
        N = 10
        res_left = self.residual[:N]
        res_right = self.residual[-N:]
       

        res_left_avg = sum(res_left)/N
       
        res_right_avg = sum(res_right)/N

        residual_average = (res_left_avg + res_right_avg)/2
       

        for j in range(2, len(self.x_array-2)):
            residual_noiseless[j] = self.residual[j] - residual_average

        max_res_location = 0
        #Finding location of top/bottom 6 maximum/minimum residual (where we want to add/remove a new peak)
        sorting_res = np.argsort(residual_noiseless)
        sorted_res = residual_noiseless[sorting_res]
       
        top_6_rev = sorted_res[-6 : ] 
        bottom_6 = sorted_res[ : 6] 
        top_6_BE = []
        bottom_6_BE = []

        top_6 = top_6_rev[::-1]
        for i in range(len(self.x_array)):
            if residual_noiseless[i] == top_6[0]:
                max_BE_1st_res = self.x_array[i]
            elif residual_noiseless[i] == top_6[1]:
                max_BE_2nd_res = self.x_array[i]
            elif residual_noiseless[i] == top_6[2]:
                max_BE_3rd_res = self.x_array[i]
            elif residual_noiseless[i] == top_6[3]:
                max_BE_4th_res = self.x_array[i]
            elif residual_noiseless[i] == top_6[4]:
                max_BE_5th_res = self.x_array[i]
            elif residual_noiseless[i] == top_6[5]:
                max_BE_6th_res = self.x_array[i]
            elif residual_noiseless[i] == bottom_6[0]:
                min_BE_1st_res = self.x_array[i]
            elif residual_noiseless[i] == bottom_6[1]:
                min_BE_2nd_res = self.x_array[i]
            elif residual_noiseless[i] == bottom_6[2]:
                min_BE_3rd_res = self.x_array[i]
            elif residual_noiseless[i] == bottom_6[3]:
                min_BE_4th_res = self.x_array[i]
            elif residual_noiseless[i] == bottom_6[4]:
                min_BE_5th_res = self.x_array[i]
            elif residual_noiseless[i] == bottom_6[5]:
                min_BE_6th_res = self.x_array[i]
        top_6_BE.append(max_BE_1st_res)
        top_6_BE.append(max_BE_2nd_res)
        top_6_BE.append(max_BE_3rd_res)
        top_6_BE.append(max_BE_4th_res)
        top_6_BE.append(max_BE_5th_res)
        top_6_BE.append(max_BE_6th_res)
        bottom_6_BE.append(min_BE_1st_res)
        bottom_6_BE.append(min_BE_2nd_res)
        bottom_6_BE.append(min_BE_3rd_res)
        bottom_6_BE.append(min_BE_4th_res)
        bottom_6_BE.append(min_BE_5th_res)
        bottom_6_BE.append(min_BE_6th_res)

        pass_cond_max = 0
        pass_cond_min = 0
        new_max_BE = 0
        new_min_BE = 0
        new_max_res = 0
        new_min_res = 0
        upper = BE_guess[0] + 2
        lower = BE_guess[0] - 2
       
        for i in range(len(top_6_BE)):
            
            if top_6_BE[i] > upper or top_6_BE[i] < lower:
                if pass_cond_max > 0:
                    pass
                else:
                    new_max_BE = top_6_BE[i]
                    new_max_res = top_6[i]
                    pass_cond_max +=1
            if bottom_6_BE[i] > upper or bottom_6_BE[i] < lower:
                if pass_cond_min > 0:
                    pass
                else:
                    new_min_BE = bottom_6_BE[i]
                    new_min_res = bottom_6[i]
                    pass_cond_min += 1


        
        

        max_peak = amp_guess[0]
        
        peak_location = 0
        for i in amp_guess: #Finding the location of the maximum amplitude peak
            if i >= max_peak:
                max_peak_location = peak_location
            peak_location += 1
        print("MINIMUM PEAK RESDIUAL", new_min_res, "-30%",-0.3*amp_guess[max_peak_location])
        print("MAXIMUM PEAK RESDIUAL", new_max_res, "15%",0.15*amp_guess[max_peak_location])
        #Check to remove peak. Won't remove if only one peak    
        if len(BE_guess) > 1:    
            if new_min_res < -0.3*amp_guess[max_peak_location] and new_min_res != 0: #NEED TO FIND A BETTER CONDITION FOR ADDING AND REMOVING PEAKS 
                #Do we want to check where the peak is located and make some evaluation from this using new_min_BE?
                print("Data is overfit. Removing peak to try to improve fitting")
                self.removePeak()

        if new_max_res > 0.15*amp_guess[max_peak_location]: #NEED TO FIND A BETTER CONDITION FOR ADDING AND REMOVING PEAKS 

            new_peak_BE = new_max_BE
            new_peak_amp = new_max_res
            print("NEW PEAK BE", new_peak_BE, "NEW PEAK AMPLITUDE", new_peak_amp)
            self.addPeak(new_peak_BE, new_peak_amp)
        else:
            print("NO NEW PEAK ADDED")
            

        
    

    

     

       


    #THIS FUNCTION IS ONLY WORKING FOR VOIGT RIGHT NOW. IMPROPER ADDITION OF THE NEW PARAMETERS WITH OTHER TYPES
    #Where do we call this function? --> Right now it is only functionable in generateFirstGen()
    def addPeak(self, new_BE, new_amp): #New function to be called after so many generations (20?). Used before next generation is created.
        """
        Appends a new peak (binding energy, amplitude) to the fit's
        parameter guesses/ranges. Part of the known-broken
        `peak_add_remove` path (crashes on the first added peak) — see
        `scanResidual`'s docstring.
        """

        BE_guess.append(new_BE) #Want to use reisdual to determine where to add new peak?
        #Condition if the new peak trying to be added is too close to another peak --> Instead adjust amp of that peak
        #How to use reisudal to find these parameters? SOS? Br? etc. --> Take on same values as other peak but give a wider range?
      
        amp_guess.append(new_amp)
        sigma_guess.append(0.5)
        gamma_guess.append(0.25)
        asym_guess.append(1)
        sos_guess.append(0)
        br_guess.append(0.5)

        amp_limited.append(False)
        BE_limited.append(False)
        gamma_limited.append(False)
        sigma_limited.append(False)
        amp_min_max = new_amp*0.2
        amp_range_min.append(-amp_min_max)
        amp_range_max.append(amp_min_max)
        amp_range_delta.append(0.05)

        BE_range_min.append(-0.2)
        BE_range_max.append(0.2)
        BE_range_delta.append(0.01)

        #Allowing for larger sigma and gamma range --> Need to figure out how to get a better estimate on these parameters from residual
        sigma_range_min.append(-0.5)
        sigma_range_max.append(0.5)
        sigma_range_delta.append(0.001)
       
        gamma_range_min.append(-0.25)
        gamma_range_max.append(0.25)
        gamma_range_delta.append(0.001)

        sos_range_min.append(-0.02)
        sos_range_max.append(0.02)
        sos_range_delta.append(0.01)

        br_range_min.append(-0.05)
        br_range_max.append(0.05)
        br_range_delta.append(0.01)

        asym_range_min.append(1)
        asym_range_max.append(1)
        asym_range_delta.append(0.01)

        
        
        peak_type.append(peak_type[0]) #Appending whatever is the 1st input peak type --> Want to change in the future to have the algorithm decide what is the best peak type
        self.peak_type = peak_type
      

        i = len(BE_guess)

        BE_correlated.append(i)
        amp_correlated.append(i)
        asym_correlated.append(i)
        sos_correlated.append(i)
        br_correlated.append(i)
        gamma_correlated.append(i)
        sigma_correlated.append(i)

        #How to allow for peak correlation for added peak? Correlate if other peaks have correlation?
        BE_correlated_mult.append(1)
        asym_correlated_mult.append(1)
        sos_correlated_mult.append(1)
        br_correlated_mult.append(1)
        amp_correlated_mult.append(1)
        sigma_correlated_mult.append(1)
        gamma_correlated_mult.append(1)

       
        #Take on the same values as the first peak? Maximum amp peak?
        #branching_ratio.append(branching_ratio[0])
        #spinOrbitSplit_guess.append(spinOrbitSplit_guess[0])
        gamma_CK_guess.append(gamma_guess[len(BE_guess)-1]) #taking in same gamma value as the newest gamma added
        asym_CK_guess.append(asym_guess[len(BE_guess)-1])
        is_singlet.append(is_singlet[0])
        is_coster_kronig.append(is_coster_kronig[0])


        npaths = len(BE_guess)
       


        #HAVE TO RECALL ALL THIS SO THAT WE CAN CALL PARS RANGE AGAIN
        photoline_select = self.photoelectronLine
        N = 10
        y_left = self.y_array[:N]
        y_right = self.y_array[-N:]
        x_left = self.x_array[:N]
        x_right = self.x_array[-N:]

        self.y_left_avg = sum(y_left)/N
        x_left_avg = sum(x_left)/N
        self.y_right_avg = sum(y_right)/N
        x_right_avg = sum(x_right)/N
        if self.y_left_avg > self.y_right_avg:
            #Should we scale baseline so that the last point (lowest BE/highest KE) is subtratced off the y values (goes to zero) and set highest possible baseline value to the left side of the data?
            baseline_range = [-(self.y_left_avg+np.round(self.y_left_avg*0.2, 2)),self.y_left_avg+np.round(self.y_left_avg*0.05, 2),1] #Changed minimum to allow for background lower. Now allows for 20% below avg. and 5 % above
           
        else:
            baseline_range = [-(self.y_left_avg+np.round(self.y_left_avg*0.2, 2)),self.y_right_avg+np.round(self.y_right_avg*0.05, 2),1]

        amp_range = []
        asym_range = []
        sos_range = []
        br_range = []
        BE_range = []
        sigma_range = []
        gamma_range = []

        max_peak = amp_range_max[0] #Finding the peak with the largest amplitude and giving it a larger BE range
        #In the future we will need to have the algorithm check the elemental peak assingment here and give it a range based on a database of collected peak information
        peak_location = 0
        for i in amp_range_max:
            if i >= max_peak:
                max_peak_location = peak_location
            peak_location += 1
        BE_range_min[max_peak_location] = -0.5
        BE_range_max[max_peak_location] = +0.5

        for i in range(len(BE_guess)): #Creating 2D arrays of range values for each peaks BE and amplitude
            
            if amp_limited[i] == True:
                amp_range.append([-0.00001, 0.00001, 0.00001]) #Making the range basically zero if the parameter is 'fixed'
            else:
                amp_range.append([amp_range_min[i], amp_range_max[i], amp_range_delta[i]])

            if asym_limited[i] == True:
                asym_range.append([-0.00001, 0.00001, 0.00001]) #Making the range basically zero if the parameter is 'fixed'
            else:
                asym_range.append([asym_range_min[i], asym_range_max[i], asym_range_delta[i]])

        



            if BE_limited[i] == True:
                
                BE_range.append([-0.00001, 0.00001, 0.00001]) #Making the range basically zero if the parameter is 'fixed'
            else:
                BE_range.append([BE_range_min[i], BE_range_max[i], BE_range_delta[i]])
        


            if sigma_limited[i] == True:
                sigma_range.append([-0.00001, 0.00001, 0.00001]) #Making the range basically zero if the parameter is 'fixed'
            else:
                sigma_range.append([sigma_range_min[i], sigma_range_max[i], sigma_range_delta[i]])



            if gamma_limited[i] == True:
                gamma_range.append([-0.00001, 0.00001, 0.00001]) #Making the range basically zero if the parameter is 'fixed'
            else:
                gamma_range.append([gamma_range_min[i], gamma_range_max[i], gamma_range_delta[i]])

            if sos_limited[i] == True:
                sos_range.append([-0.00001, 0.00001, 0.00001]) #Making the range basically zero if the parameter is 'fixed'
            else:
                sos_range.append([sos_range_min[i], sos_range_max[i], sos_range_delta[i]])

            if br_limited[i] == True:
                br_range.append([-0.00001, 0.00001, 0.00001]) #Making the range basically zero if the parameter is 'fixed'
            else:
                br_range.append([br_range_min[i], br_range_max[i], br_range_delta[i]])


        #REDFINE parse_range after appending everything needed to add a peak


        #Recalling self.pars_range after fitting has already begun seems to restart the guesses for values to the values the algorithm has found and creates new ranges around those values
        self.pars_range = {
            'Binding Energy': BE_range,
            'BE':BE_guess,
            'BE_limited' :BE_limited,
            'BE_correlated' :BE_correlated,
            'BE_correlated_mult' :BE_correlated_mult,
            'Gaussian': sigma_range,
            'Sigma': sigma_guess,
            'sigma_limited' :sigma_limited,
            'sigma_correlated' :sigma_correlated,
            'sigma_correlated_mult' :sigma_correlated_mult,
            'Lorentzian': gamma_range,
            'Gamma': gamma_guess,
            'gamma_limited' :gamma_limited,
            'gamma_correlated' :gamma_correlated,
            'gamma_correlated_mult' :gamma_correlated_mult,
            'Amplitude' : amp_range,
            'Amp': amp_guess,
            'amp_limited' :amp_limited,
            'amp_correlated' :amp_correlated,
            'amp_correlated_mult' :amp_correlated_mult,
            'Asymmetry' : asym_range,
            'Asym': asym_guess,
            'asym_limited' :asym_limited,
            'asym_correlated' : asym_correlated,
            'asym_correlated_mult' :asym_correlated_mult,
            #'Asymmetry' : asymmetry_range,
            'Asymmetry Doniach-Sunjic' : asymmetryDoniach_range,
            'k_range' : k_range,
            'Background' : background_range,
            #'Shirley Background' : background_shir_range,
            'Slope' : slope_range,
            'npeaks' : npaths,
            'baseline' : baseline_range,
            #'branching_ratio_range' : branching_ratio_range,
            #'branching_ratios' : branching_ratio,
            'br_range' : br_range,
            'BR': br_guess,
            'br_limited' :br_limited,
            'br_correlated' : br_correlated,
            'br_correlated_mult' :br_correlated_mult,
            'is_singlet' : is_singlet,
            'is_coster_kronig' : is_coster_kronig,
            'Lorentzian Coster-Kronig' : gamma_CK_range,
            'Gamma Coster-Kronig' : gamma_CK_guess,
            'Asymmetry Range Coster-Kronig' : asym_CK_range,
            'Asymmetry Coster-Kronig' : asym_CK_guess,
            #'spinOrbitSplitting' : spinOrbitSplit_range,
            #'spinOrbitSplit' : spinOrbitSplit_guess, #Are we still using this one??
            'spinOrbitSplitting' : sos_range,
            'SOS': sos_guess,
            'sos_limited' :sos_limited,
            'sos_correlated' : sos_correlated,
            'sos_correlated_mult' :sos_correlated_mult,
            'photoline' : photoline_select
        }

    def removePeak(self): #Removing peak from array. Right now it just removes the last input --> Need to change that
        """
        Drops the last peak's guesses/ranges (not necessarily the worst
        one — see the TODO in the signature comment). Part of the
        known-broken `peak_add_remove` path — see `scanResidual`'s
        docstring.
        """

        BE_guess.pop()
      
        amp_guess.pop()
        sos_guess.pop()
        br_guess.pop()
        asym_guess.pop()
        sigma_guess.pop()
        gamma_guess.pop()

        amp_limited.pop()
        asym_limited.pop()
        sos_limited.pop()
        br_limited.pop()
        BE_limited.pop()
        gamma_limited.pop()
        sigma_limited.pop()
        
        amp_range_min.pop()
        amp_range_max.pop()
        amp_range_delta.pop()

        sos_range_min.pop()
        sos_range_max.pop()
        sos_range_delta.pop()

        br_range_min.pop()
        br_range_max.pop()
        br_range_delta.pop()

        asym_range_min.pop()
        asym_range_max.pop()
        asym_range_delta.pop()

        BE_range_min.pop()
        BE_range_max.pop()
        BE_range_delta.pop()

       
        sigma_range_min.pop()
        sigma_range_max.pop()
        sigma_range_delta.pop()

        gamma_range_min.pop()
        gamma_range_max.pop()
        gamma_range_delta.pop()

        peak_type.pop()
        self.peak_type = peak_type
      

       

        BE_correlated.pop()
        amp_correlated.pop()
        asym_correlated.pop()
        sos_correlated.pop()
        br_correlated.pop()
        gamma_correlated.pop()
        sigma_correlated.pop()
       
        BE_correlated_mult.pop()
        amp_correlated_mult.pop()
        asym_correlated_mult.pop()
        sos_correlated_mult.pop()
        br_correlated_mult.pop()
        sigma_correlated_mult.pop()
        gamma_correlated_mult.pop()


        #branching_ratio.pop()
        #spinOrbitSplit_guess.pop()
        gamma_CK_guess.pop()
        asym_CK_guess.pop()
        is_singlet.pop()
        is_coster_kronig.pop()
        npaths = len(BE_guess)
       


        #HAVE TO RECALL ALL THIS SO THAT WE CAN CALL PARS RANGE AGAIN
        photoline_select = self.photoelectronLine
        N = 10
        y_left = self.y_array[:N]
        y_right = self.y_array[-N:]
        x_left = self.x_array[:N]
        x_right = self.x_array[-N:]

        self.y_left_avg = sum(y_left)/N
        x_left_avg = sum(x_left)/N
        self.y_right_avg = sum(y_right)/N
        x_right_avg = sum(x_right)/N
        if self.y_left_avg > self.y_right_avg:
            #Should we scale baseline so that the last point (lowest BE/highest KE) is subtratced off the y values (goes to zero) and set highest possible baseline value to the left side of the data?
            baseline_range = [-(self.y_left_avg+np.round(self.y_left_avg*0.2, 2)),self.y_left_avg+np.round(self.y_left_avg*0.05, 2),1] #Changed minimum to allow for background lower. Now allows for 20% below avg. and 5 % above
           
        else:
            baseline_range = [-(self.y_left_avg+np.round(self.y_left_avg*0.2, 2)),self.y_right_avg+np.round(self.y_right_avg*0.05, 2),1]

        amp_range = []
        sos_range = []
        br_range = []
        asym_range = []
        BE_range = []
        sigma_range = []
        gamma_range = []

        max_peak = amp_range_max[0] #Finding the peak with the largest amplitude and giving it a larger BE range
        #In the future we will need to have the algorithm check the elemental peak assingment here and give it a range based on a database of collected peak information
        '''
        peak_location = 0
        for i in amp_range_max:
            if i >= max_peak:
                max_peak_location = peak_location
            peak_location += 1
        BE_range_min[max_peak_location] = -0.5
        BE_range_max[max_peak_location] = +0.5
        '''

        for i in range(len(BE_guess)): #Creating 2D arrays of range values for each peaks BE and amplitude
          
            if amp_limited[i] == True:
                amp_range.append([-0.00001, 0.00001, 0.00001]) #Making the range basically zero if the parameter is 'fixed'
            else:
                amp_range.append([amp_range_min[i], amp_range_max[i], amp_range_delta[i]])

            if asym_limited[i] == True:
                asym_range.append([-0.00001, 0.00001, 0.00001]) #Making the range basically zero if the parameter is 'fixed'
            else:
                asym_range.append([asym_range_min[i], asym_range_max[i], asym_range_delta[i]])



            if BE_limited[i] == True:
                
                BE_range.append([-0.00001, 0.00001, 0.00001]) #Making the range basically zero if the parameter is 'fixed'
            else:
                BE_range.append([BE_range_min[i], BE_range_max[i], BE_range_delta[i]])
        


            if sigma_limited[i] == True:
                sigma_range.append([-0.00001, 0.00001, 0.00001]) #Making the range basically zero if the parameter is 'fixed'
            else:
                sigma_range.append([sigma_range_min[i], sigma_range_max[i], sigma_range_delta[i]])



            if gamma_limited[i] == True:
                gamma_range.append([-0.00001, 0.00001, 0.00001]) #Making the range basically zero if the parameter is 'fixed'
            else:
                gamma_range.append([gamma_range_min[i], gamma_range_max[i], gamma_range_delta[i]])

            if sos_limited[i] == True:
                sos_range.append([-0.00001, 0.00001, 0.00001]) #Making the range basically zero if the parameter is 'fixed'
            else:
                sos_range.append([sos_range_min[i], sos_range_max[i], sos_range_delta[i]])

            if br_limited[i] == True:
                br_range.append([-0.00001, 0.00001, 0.00001]) #Making the range basically zero if the parameter is 'fixed'
            else:
                br_range.append([br_range_min[i], br_range_max[i], br_range_delta[i]])


        #REDFINE parse_range after appending everything needed to add a peak
      

        self.pars_range = {
            'Binding Energy': BE_range,
            'BE':BE_guess,
            'BE_limited' :BE_limited,
            'BE_correlated' :BE_correlated,
            'BE_correlated_mult' :BE_correlated_mult,
            'Gaussian': sigma_range,
            'Sigma': sigma_guess,
            'sigma_limited' :sigma_limited,
            'sigma_correlated' :sigma_correlated,
            'sigma_correlated_mult' :sigma_correlated_mult,
            'Lorentzian': gamma_range,
            'Gamma': gamma_guess,
            'gamma_limited' :gamma_limited,
            'gamma_correlated' :gamma_correlated,
            'gamma_correlated_mult' :gamma_correlated_mult,
            'Amplitude' : amp_range,
            'Amp': amp_guess,
            'amp_limited' :amp_limited,
            'amp_correlated' :amp_correlated,
            'amp_correlated_mult' :amp_correlated_mult,
            'Asymmetry' : asym_range,
            'Asym': asym_guess,
            'asym_limited' :asym_limited,
            'asym_correlated' : asym_correlated,
            'asym_correlated_mult' :asym_correlated_mult,
            #'Asymmetry' : asymmetry_range,
            'Asymmetry Doniach-Sunjic' : asymmetryDoniach_range,
            'k_range' : k_range,
            'Background' : background_range,
            #'Shirley Background' : background_shir_range,
            'Slope' : slope_range,
            'npeaks' : npaths,
            'baseline' : baseline_range,
            #'branching_ratio_range' : branching_ratio_range,
            #'branching_ratios' : branching_ratio,
            'br_range' : br_range,
            'BR': br_guess,
            'br_limited' :br_limited,
            'br_correlated' : br_correlated,
            'br_correlated_mult' :br_correlated_mult,
            'is_singlet' : is_singlet,
            'is_coster_kronig' : is_coster_kronig,
            'Lorentzian Coster-Kronig' : gamma_CK_range,
            'Gamma Coster-Kronig' : gamma_CK_guess,
            'Asymmetry Range Coster-Kronig' : asym_CK_range,
            'Asymmetry Coster-Kronig' : asym_CK_guess,
            #'spinOrbitSplitting' : spinOrbitSplit_range,
            #'spinOrbitSplit' : spinOrbitSplit_guess, #Are we still using this one??
            'spinOrbitSplitting' : sos_range,
            'SOS': sos_guess,
            'sos_limited' :sos_limited,
            'sos_correlated' : sos_correlated,
            'sos_correlated_mult' :sos_correlated_mult,
            'photoline' : photoline_select
        }
       
        
        




    def generateIndividual(self):
        """
        Generate singular individual
        """
        
        ind = Individual(self.backgrounds,peak_type,self.scale_var,self.pars_range)
        return ind

    def generateFirstGen(self):
        print("generating first gen")
        self.Populations=[]

        #ind = Individual(self.backgrounds,peak_type,self.scale_var,self.pars_range)
        #self.addPeak(ind)

        
        for i in range(self.npops):
          
            self.Populations.append(self.generateIndividual())
        self.eval_Population()
        self.globBestFit = self.sorted_population[0]

        print("First gen generated")
    # @profile
    def fitness(self,indObj):
        """
        Evaluate fitness of an individual

        The loss math lives in xps_neo.loss.compute_loss, shared with the
        GUI post-analysis (Phase 3 unification).
        """
        Individual = indObj

        #HOW TO MAKE SCALE VALUE SEEN HERE?
        yTotal = Individual.getFit(self.x_array,self.y_array, self.backgrounds)

        loss, self.residual = compute_loss(
            self.x_array, self.y_array, yTotal,
            self.y_left_avg, self.y_right_avg,
            apply_penalties=self.genNum < self.gen_alt)
        return loss, yTotal

    def _eval_many(self, populations):
        """Evaluate a list of individuals, serial or parallel (--workers).

        Fitness evaluation consumes no RNG (all draws happen in the
        mutation/crossover steps in this process), so parallel evaluation
        follows the identical trajectory. Workers return the individuals
        because getFit writes evaluated values back onto peak attributes;
        the caller must replace its references with the returned objects.

        Returns a list of (individual, loss, yTotal, residual) in
        population order.
        """
        if getattr(self, "workers", 1) > 1:
            if getattr(self, "_pool", None) is None:
                import multiprocessing
                self._pool = multiprocessing.Pool(self.workers)
            apply_pen = self.genNum < self.gen_alt
            jobs = [(ind, self.x_array, self.y_array, self.y_left_avg,
                     self.y_right_avg, apply_pen, self.backgrounds)
                    for ind in populations]
            results = []
            for orig, (evaluated, loss, yTotal, residual) in zip(
                    populations, self._pool.map(_evaluate_individual, jobs)):
                _absorb_state(orig, evaluated)
                results.append((orig, loss, yTotal, residual))
            return results
        results = []
        for ind in populations:
            loss, yTotal = self.fitness(ind)
            results.append((ind, loss, yTotal, self.residual))
        return results
   
    def eval_Population(self):
        """
        Evalulate populations
        """

        score = []
        populationPerf = {}
        self.og_fitness = []
        self.nan_counter = 0
        best_y = []
        temp_temp_score = 10000000000000000000000000
        yPerf = {}
        self.nan_counter = 0
        
        for i, (individual, temp_score, bestY, residual) in enumerate(
                self._eval_many(self.Populations)):
            self.Populations[i] = individual  # returned copy carries getFit's attribute writes
            self.residual = residual

            #Finding the best y values from the best individual in the fit so far
            if temp_score < temp_temp_score:
                best_y = bestY
                temp_temp_score = temp_score
            # Calculate the score, if encounter nan, discard and generate new individual later
            if np.isnan(temp_score): #Change to be if residual is bad AND numGenSinceImproved > 10% of generation size then .... add/remove peak OR adjust parameters
                self.nan_counter +=1
            else:
                score.append(temp_score)
                populationPerf[individual] = temp_score

        #Added this here to redefine param ranges into a 2D array --> Needs to be 1D to be read in xps_fit.py
        '''
        self.full_amp_range = self.pars_range['Amplitude']
        self.full_asym_range = self.pars_range['Asymmetry']
        self.full_sigma_range = self.pars_range['Gaussian']
        self.full_gamma_range = self.pars_range['Lorentzian']
        self.full_energy_range = self.pars_range['Binding Energy']
        
       
        #Individual.defineRanges(self)
        '''
                
                
        self.og_fitness = score
        self.sorted_population = sorted(populationPerf.items(), key=operator.itemgetter(1), reverse=False) #Get an error here if y values go below zero for some reason
     
       
        '''
        for a,b in self.sorted_population:
            print(str(b) + " " + str(a.get_params()))
        '''
        ''' Debugging again
        for i in range(len(self.sorted_population)):
            print(self.sorted_population[i][0].get_peak(0))
        '''
        self.currBestFit = self.sorted_population[0]
        self.currBestY = best_y
      

        return score
    
    def eval_Pop(self,populations):
        scores = []
        self.nan_counter = 0
        for i, (individual, temp_score, _yTotal, residual) in enumerate(
                self._eval_many(populations)):
            populations[i] = individual  # returned copy carries getFit's attribute writes
            self.residual = residual
            if np.isnan(temp_score):
                self.nan_counter +=1
            else:
                scores.append(temp_score)

            #scores.append(temp_score)

        return scores

    
    def next_generation(self):
        """
        Calculate next generations

        """
        mut_De_counter = 0
        self.st = time.time()
        # ray.init()
        self.logger.info("---------------------------------------------------------")
        self.logger.info(datetime.datetime.fromtimestamp(self.st).strftime('%Y-%m-%d %H:%M:%S'))
        self.logger.info(f"{bcolors.BOLD}Gen: {bcolors.ENDC}{self.genNum+1}")

        self.genNum += 1

        # Evaluate Fittness
        score = self.eval_Population()
        self.bestDiff = abs(self.globBestFit[1]-self.currBestFit[1])
        if self.currBestFit[1] < self.globBestFit[1]:
            self.globBestFit = self.currBestFit
            self.numGenSinceImproved = 0
        else:
            self.numGenSinceImproved += 1

        
        self.mutatePopulation()
        if self.mut_opt != 3:
            self.selectFromPopulation()
            self.createChildren()
            self.logger.info("Number of Breeders: " + str(len(self.parents)))
            self.logger.info("DiffCounter: " + str(self.diffCounter))
            self.logger.info("Diff %: " + str(self.diffCounter / self.genNum))
            self.logger.info("Mutation Chance: " + str(self.mut_chance))
        # DE
        else:
            self.crossoverPopulation()
            self.adjust_DE_parameters()
            trial_fitness = self.eval_Pop(self.trialPopulations)
           
            #print(trial_fitness[0], self.og_fitness[0])
            for i in range(self.npops): #Added in try and except because error would not go away
                try:
                    if trial_fitness[i] < self.og_fitness[i]: #This line keeps raising an error for when it is not true --> Try hitting Run again and it should work. If not close and reopen gui.

                        self.Populations[i] = self.trialPopulations[i]
                        mut_De_counter += 1
                    else:
                        pass
                except:
                    pass
            self.logger.info("Mutated Individuals Added to Population: " + str(mut_De_counter))    

            self.logger.info(f"Average Trial Population Fitness: {np.average(trial_fitness)}") #Don't really think this values matters since trial population is all over the place
            self.logger.info(f"Average Population Fitness: {np.average(self.og_fitness)}")

        self.et = timecall()
        self.tdiff = self.et - self.st
        self.tt = self.tt + self.tdiff

        self.verboseGeneration()
        self.logger.info("Time: "+ str(round(self.tdiff,5))+ "s")


    def verboseGeneration(self):
        """Print the best fit in a generation and
        """
        with np.printoptions(precision=5, suppress=True):
            self.logger.info("Different from last best fit: " +str(self.bestDiff))
            self.logger.info("Num of Generation since improved: " + str(self.numGenSinceImproved))
            self.logger.info(bcolors.BOLD + "Best fit: " + bcolors.OKBLUE + str(self.currBestFit[1]) + bcolors.ENDC)
            self.logger.info("Best fit combination:\n" + str((self.sorted_population[0][0].get_params()))) #This is being seen to have doublet parameters but xps_analysis is not able to see them...unsure why
            self.logger.info(bcolors.BOLD + "History Best: " + bcolors.OKBLUE + str(self.globBestFit[1]) + bcolors.ENDC)
            # self.logger.info("NanCounter: " + str(self.nan_counter))
            self.logger.info("History Best Indi:\n" + str((self.globBestFit[0].get_params())))

      


    def crossoverPopulation(self):
        """DE crossover pass: builds `trialPopulations` by crossing each
        mutated individual with its corresponding population member
        (`crossoverDE`), used when `mutated_options == 3` (DE mode)."""
        self.trialPopulations = []
        for i in range(self.npops):
            self.trialPopulations.append(self.crossoverDE(self.mutated_Populations[i],self.Populations[i],self.cR))


    def crossoverDE(self,mutateInd: Individual, popInd: Individual, cR: float) -> Individual:
        """Uniform crossover of DE using `self.cR`

        Args:
            mutateInd (Individual): Individual for the mutated population
            popInd (Individual): Individual for the original population
            cR (float): crossover rate

        Returns:
            Individual: crossovered individual
        """

        p = np.random.rand(len(mutateInd))
        mutatePars = mutateInd.get_params()
        popPars = popInd.get_params()
        tempPars = []
        for i in range(len(mutateInd)):
            if p[i] < cR:
                tempPars.append(mutatePars[i])
            else:
                tempPars.append(popPars[i])

        split_full_list = XPS_GA.split_into_x(tempPars)

        temp_individual = self.generateIndividual()
        XPS_GA.setPars(temp_individual,split_full_list)

        return temp_individual

    def adjust_DE_parameters(self,on=True):
        """
        Adjust the DE parameters using jDE algorithm

        Args:
            on (bool, optional): _description_. Defaults to True.
        """
        if on:
            rand_val = np.random.rand(4)
            tau_1 = 0.1
            tau_2 = 0.1
            if rand_val[1] < tau_1:
                self.F = 0.1 + rand_val[0] * 0.9
                self.logger.info(f"F has been adjusted to {np.round(self.F,4)}")

            if rand_val[3] < tau_2:
                self.cR = rand_val[2]
                self.logger.info(f"Cr has been adjusted to {np.round(self.cR,4)}")
       
    def mutatePopulation(self):
        """
        # Mutation operators
        # 0 = original: generated a new versions: #Rechenberg Mutation 
        # 1 = mutated every genes in the total populations #Random Perturbations
        # 2 = mutated genes inside population based on secondary probability #Metropolis Mutation 
        # 3 = Differential Evolution

        # TODO:
            options 2 and 3 needs to reimplmented
        """
        self.nmutate = 0
    
        
        "DE mutation"
        if self.mut_opt == 3:
            self.mutated_Populations = []
            for i in range(self.npops):
                candidates = [candidate for candidate in range(self.npops) if candidate != i]
                a,b,c = np.random.choice(candidates,3,replace=False)
                mutation_vectors = [self.Populations[a],self.Populations[b],self.Populations[c]]
                temp_individual = self.mutateIndi_DE(mutation_vectors,self.F)
                temp_individual.checkBound()
                self.mutated_Populations.append(temp_individual)

        else: 
            # Rechenberg mutation
            if self.mut_opt  == 0:
                if self.genNum > 20:
                    if self.bestDiff < 0.1:
                        self.diffCounter += 1
                    else:
                        self.diffCounter -= 1
                    if (abs(self.diffCounter)/ float(self.genNum)) > 0.2:
                        self.mut_chance += 0.5
                        self.mut_chance = abs(self.mut_chance)
                    elif (abs(self.diffCounter) / float(self.genNum)) < 0.2:
                        self.mut_chance -= 0.5
                        self.mut_chance = abs(self.mut_chance)

            for i in range(self.npops):
                if random.random()*100 < self.mut_chance:
                    self.nmutate += 1
                    self.Populations[i] = self.mutateIndi(i)

            self.logger.info("Mutate Times: " + str(self.nmutate)) #Moved inside else becuase this doesnt change for DE mutation option. Meaningless 

    def setPars(individual,split_full_list):
        """Static Method to set the parameters

        Args:
            individual (individuals): individuals type
            split_full_list (list): multidims list containing each peaks type and
                parameters, and bg peaks.
        """
        full_list = copy.copy(split_full_list)
        # remove the bg
        bg = full_list.pop()
        for i in range(len(full_list)):
            individual.setPeak(i,full_list[i])

    @staticmethod
    def split_into_x(full_list):
        """Split a full parameters list into multidmenisonal list

        Args:
            full_list (list): full list of parameters
        """
        split_list = []
        temp_list = []
        for i in full_list:
            if isinstance(i,str) == False:
                temp_list.append(i)
            else:
                temp_list.append(i)
                split_list.append(temp_list)
                temp_list = []

        return split_list

    def mutateIndi_DE(self,mutateIndividuals:list,F:float) -> Individual:
        """
        Mutate the individuals using DE mutation

        Args:
            mutated_individuals (list): Input list of individuals
            F (float): Mutation Factors
        """
        length = len(mutateIndividuals[0])
        assert all(len(lst) == length for lst in mutateIndividuals)

        x_Pars = mutateIndividuals[0].get_params()
        y_Pars = mutateIndividuals[1].get_params()
        z_Pars = mutateIndividuals[2].get_params()

        # Convert to the extracted full list with peaks and BG
        full_list = []
        for i in range(len(x_Pars)):
           
            if isinstance(x_Pars[i],str) == False:
                full_list.append(x_Pars[i]  + F*(y_Pars[i]  - z_Pars[i]))
            else:
                full_list.append(x_Pars[i])

        split_full_list = XPS_GA.split_into_x(full_list)
        temp_individual = self.generateIndividual()
   
        XPS_GA.setPars(temp_individual,split_full_list)

        return temp_individual

    def mutateIndi(self,indi):
        """
        Generate new individual during mutation operator
        """
        if self.mut_opt == 0:
            # Create a new individual with Rechenberg
            newIndi = self.generateIndividual()
        # Random pertubutions
        if self.mut_opt == 1:
            # Random Pertubutions
            self.Populations[indi].mutate_(self.mut_chance)
            newIndi = self.Populations[indi]
            # Mutate every gene in the Individuals

        if self.mut_opt == 2:
            #Metropoli mutation 
            n_success = 0
            og_individual = self.generateIndividual()
            # Create a new individual with the same parameters
            og_pars = copy.copy(self.Populations[indi].get_peaks())

            #og_individual.setPeaks(og_pars)
            og_score, yTot = self.fitness(og_individual)

            new_individual = self.generateIndividual()
            mut_score, yTot_mut = self.fitness(new_individual)

            T = - self.bestDiff/(np.log(1-(self.genNum/self.ngen)+self.tol))
           
            if mut_score < og_score:
                n_success = n_success + 1

                newIndi = new_individual
            elif np.exp(-(mut_score-og_score)/(T+self.tol)) > np.random.uniform():
                n_success = n_success + 1
                newIndi = new_individual
            else:
                newIndi = og_individual

        return newIndi
        '''
        if self.mut_opt == 3: #Not working --> Replace with DE?
            def delta_fun(t,delta_val):
                rnd = np.random.random()
                return delta_val*(1-rnd**(1-(t/self.ngen))**5)

            og_indi = copy.deepcopy(self.Populations[indi])
            og_data = og_indi.get_var() #error --> no value get_var
            for i,path in enumerate(og_data):
                print(i,path)
                arr = np.random.randint(2,size=3)
                for j in range(len(arr)):
                    new_path = []
                    val = path[j]
                    if arr[j] == 0:
                        UP = self.pathrange_Dict[i].get_lim()[j+1][1]
                        del_val = delta_fun(self.genNum,UP-val)
                        val = val + del_val
                    if arr[j] == 1:
                        LB = self.pathrange_Dict[i].get_lim()[j+1][0]
                        del_val = delta_fun(self.genNum,val-LB)
                    new_path.append(val)
                self.Populations[indi].set_path(i,new_path[0],new_path[1],new_path[2])
        if self.mut_opt == 4:
            newIndi = self.generateIndividual(self.bestE0)
        '''
        return newIndi

    def selectFromPopulation(self):
        self.parents = []

        select_val = np.minimum(self.n_bestsam,len(self.sorted_population))
        self.n_recover = 0
        if len(self.sorted_population) < self.n_bestsam:
            self.n_recover = self.n_bestsam - len(self.sorted_population)
        for i in range(select_val):
            self.parents.append(self.sorted_population[i][0])

    @staticmethod
    def crossoverParams(val1,val2):
        if np.random.randint(0,2):
            return val1
        else:
            return val2

    def crossover(self,individual1, individual2):
        """
        Uniform Cross-Over, 50% percentage chance
        """
        child = self.generateIndividual()

        individual1_path = individual1.get_params()
        individual2_path = individual2.get_params()

        #print("Ind 1 : " + str(individual1_path))
        #print("Ind 2 : " + str(individual2_path))
        temp_path = []
        dividers = [] # markers where the strings are in the list of params, this indicates where the array switches to a new peak or background
        #crossover for peak vars
        for j in range(len(individual1_path)):
            if (isinstance(individual1_path[j],str)):
                dividers.append(j)
            if np.random.randint(0,2) == True:
                temp_path.append(individual1_path[j])
            else:
                temp_path.append(individual2_path[j])
            '''
        for j in range(1):
            if np.random.randint(0,2) == True:
                temp_path.append(individual1_path[1][j])
            else:
                temp_path.append(individual2_path[1][j])
        '''
        #print("Temp Path: " + str(temp_path))
        temp_peak = []
        #print(temp_path)
        divider = 0
        peakNum = 0
        bkgnNum = 0
        for k in range(len(dividers)):
            for j in range(divider,dividers[k]+1):
                temp_peak.append(temp_path[j])
            if i < self.npaths:
                #print()
                #print("Child pre-write: " + str(child.get_params()))
                #print("temp peak : " + str(temp_peak))
                if child.setPeak(peakNum,temp_peak) == -1:
                    if bkgnNum<len(background_type):
                        #print("Bkgn")
                        child.setBkgn(bkgnNum,temp_peak)
                        bkgnNum += 1
                else:
                    #print("wrote peak")
                    peakNum +=1
                #print("Child after write")
                #print(child.get_params())
                #print()
                temp_peak = []
            divider = j + 1

        #print("Child : " + str(child.get_params()))
        '''
        child.setPeak(i,temp_path[0],temp_path[1],temp_path[2],temp_path[3])
        child.get_background(0).set_k(temp_path[4])
        '''
        '''
        print(temp_path)
        print("Child:")
        print(child.get_params())
        sys.exit(1)
        '''
        return child

    def createChildren(self):
        """
        Generate Children
        """
        self.nextPopulation = []
        # --- append the breeder ---
        for i in range(len(self.parents)):
            self.nextPopulation.append(self.parents[i])
        # print(len(self.nextPopulation))
        # --- use the breeder to crossover
        # print(abs(self.npops-self.n_bestsam)-self.n_lucksam)

        for i in range(abs(self.npops-self.n_bestsam)-self.n_lucksam):
            par_ind = np.random.choice(len(self.parents),size=2,replace=False)
            child = self.crossover(self.parents[par_ind[0]],self.parents[par_ind[1]])
            self.nextPopulation.append(child)
        # print(len(self.nextPopulation))

        for i in range(self.n_lucksam):
            self.nextPopulation.append(self.generateIndividual())
        # print(len(self.nextPopulation))

        for i in range(self.n_recover):
            self.nextPopulation.append(self.generateIndividual())

        # for i in range(self.nan_counter):
        #     self.nextPopulation.append(self.generateIndividual())

        random.shuffle(self.nextPopulation)
        self.Populations = self.nextPopulation

    def run_verbose_start(self):
        self.logger.info("-----------Inputs File Stats---------------")
        self.logger.info(f"{bcolors.BOLD}File{bcolors.ENDC}: {self.data_file}")
        #self.logger.info(f"{bcolors.BOLD}File Type{bcolors.ENDC}: {self.data_obj._ftype}")
        self.logger.info(f"{bcolors.BOLD}File{bcolors.ENDC}: {self.output_path}")
        self.logger.info(f"{bcolors.BOLD}Population{bcolors.ENDC}: {self.npops}")
        self.logger.info(f"{bcolors.BOLD}Num Gen{bcolors.ENDC}: {self.ngen}")
        self.logger.info(f"{bcolors.BOLD}Mutation Opt{bcolors.ENDC}: {self.mut_opt}")
        self.logger.info("-------------------------------------------")

    def run_verbose_end(self):
        self.logger.info("-----------Output Stats---------------")
        # self.logger.info(f"{bcolors.BOLD}Total)
        self.logger.info(f"{bcolors.BOLD}Total Time(s){bcolors.ENDC}: {round(self.tt,4)}")
        self.logger.info("-------------------------------------------")

    def run(self, data_peak_add):
        """Drives the full generation loop: evaluate -> select/crossover/
        mutate (GA) or crossover/mutate/select (DE, `mutated_options == 3`)
        -> re-evaluate -> log/write outputs, optionally scanning for
        peak add/remove each generation when `data_peak_add` is set."""
        self.run_verbose_start()
        self.historic = []
        self.historic.append(self.Populations)
        count = 0
        peak_add = 0
        before_best_fit = 0
        num_peaks = len(BE_guess)
        for i in range(self.ngen):
           
            temp_gen = self.next_generation()
        
       
            #Peak addition/removal determined by user checking selection button in GUI
            #if alt_lorentz == True: # Add in check to see right side of peaks (low BE/high KE). If the fit is low on the right side then increase LFWHM value. If no change THEN go to peak addition if selected as an option. 
                
            
            if peak_add_remove == True:
                
                if self.genNum >= self.gen_alt:
                #if count == 5*len(BE_guess): #Condition for adding new peak. Does so after so many generations --> Maybe make peak num dependent 10*numPeaks --> Need more time to analyze more peaks
                
                    if peak_add > 0: 
                    
                        if before_best_fit < self.currBestFit[1]: #Check to make sure that the added peak improved the fit
                            #REMOVE PEAK THAT WAS ADDED
                            self.removePeak()
                            print("No improvement from peak addition. Added peak removed")
                    
                
                    self.scanResidual()
                    if len(BE_guess) > num_peaks:
                        num_peaks = len(BE_guess)
                        self.generateFirstGen()
                        if os.path.exists(self.file_data): os.remove(self.file_data)
                    elif len(BE_guess) < num_peaks:
                        num_peaks = len(BE_guess)
                        self.generateFirstGen()
                        if os.path.exists(self.file_data): os.remove(self.file_data) #CAN ADD NEW PEAKS HERE BUT IT WILL NEED TO DESTROY THE OLD OUTPUT FILE AND CREATE A NEW ONE SO THAT THE NUMBER OF COLUMNS ARE ALL THE SAME FOR ANALYSIS
                    before_best_fit = self.currBestFit[1]
                    count = 0
                    peak_add += 1
            
            
            self.output_generations()
            count += 1
           

        #print(self.globBestFit[0].getFit(self.x_array,self.y_array))
        self.run_verbose_end()
        # test_y = self.export_paths(self.globBestFit[0])
        # plt.plot(self.data_obj.get_raw_data()[:,0],self.data_obj.get_raw_data()[:,1],'b-.')
        # plt.plot(self.x_slice,self.y_slice,'o--',label='data')
        # plt.plot(self.x_slice,test_y,'r--',label='model')
        # plt.legend()
        # plt.show()

    def export_paths(self,indObj):
        area_list=[]
        Individual = indObj.get()

        yTotal = np.zeros(len(self.x_slice))
        plt.figure()
        for i,paths in enumerate(Individual):
            y = paths.getY()

            yTotal += y
            # area = np.trapz(y.flatten(),x=self.x_slice.flatten())
            # component = paths.get_func(self.x_slice).reshape(-1,1)

            # area_list.append(area)

        Total_area = np.sum(area_list)
        return yTotal

    def output_generations(self):
        """
        Output generations result into two files
        """
        try:
            f1 = open(self.file,"a")
            f1.write(str(self.genNum) + "," + str(self.tdiff) + "," +
                str(self.currBestFit[1]) + "," + str(self.currBestFit[0].get_params()) +")," +
                str(self.globBestFit[1]) + "," + str(self.globBestFit[0].get_params()) +"\n")
        finally:
            f1.close()
        try:
            f2 = open(self.file_data,"a")
            write = csv.writer(f2)
            bestFit = self.globBestFit[0]
            #write.writerow((bestFit[i][0], bestFit[i][1], bestFit[i][2]))
            str_pars = bestFit.get_params(for_output_file = True)
            write.writerow(str_pars)
            f2.write("#################################\n")
        finally:
            f2.close()

    def __init__(self,scale_var = False,x_offset = 0.0, y_offset = 0.0, data_KE = False, data_XES = False, data_peak_add = False, workers = 1, verbose = False):
        """
        Steps to Initalize XPS
            XPS
        """
        self.workers = workers
        self.verbose = verbose
        self._pool = None
        self.scale_var = scale_var
        # initialize params
        self.initialize_params()
        # variables
        self.initialize_variable()
        # initialze file paths
        self.initialize_file_path()
        # initialize range
        self.initialize_range(scale_var,x_offset,y_offset, data_KE, data_XES)
        # Generate first generation
        self.generateFirstGen()

def main(argv=None):
    """CLI entry point (`xps_neo`): parse args -> load/validate the .ini ->
    build an XPS_GA and run it."""
    args = parse_args(argv)
    try:
        file_dict = load_file_dict(args)
        config = load_config(file_dict)
    except KeyError as exc:
        print(f"xps_neo: invalid input file: {exc.args[0]} "
              f"(see tests/golden/cases/ for a working reference)",
              file=sys.stderr)
        sys.exit(2)
    except FileNotFoundError as exc:
        print(f"xps_neo: {exc}", file=sys.stderr)
        sys.exit(2)
    # The methods of XPS_GA read config values as bare module globals (a
    # Phase 2 transitional seam; full threading comes with unification).
    # Installing the loaded config here preserves the exact pre-refactor
    # semantics: reads see these bindings, in-place mutation is shared,
    # and local rebinding inside methods stays local.
    globals().update(config)
    try:
        ga = XPS_GA(workers=args.workers, verbose=args.verbose)
    except FileNotFoundError as exc:
        print(f"xps_neo: data file not found: {exc}", file=sys.stderr)
        sys.exit(2)
    try:
        ga.run(data_peak_add=False)
    finally:
        if getattr(ga, "_pool", None) is not None:
            ga._pool.terminate()


def _evaluate_individual(args):
    """Multiprocessing worker for --workers population evaluation.

    Runs in a child process: evaluates one individual's curve and loss.
    No RNG is consumed here (verified: the shape/background math is
    deterministic), so parallel evaluation follows the GA's serial
    trajectory exactly; the individual is returned because getFit writes
    evaluated values back onto peak attributes.
    """
    (individual, x_array, y_array, y_left_avg, y_right_avg,
     apply_penalties, backgrounds) = args
    yTotal = individual.getFit(x_array, y_array, backgrounds)
    loss, residual = compute_loss(x_array, y_array, yTotal,
                                  y_left_avg, y_right_avg, apply_penalties)
    return individual, loss, yTotal, residual


# Structural/shared attributes that must keep their main-process identity
# when absorbing worker state (aliasing between Populations, globBestFit,
# parents and duplicate slots is load-bearing in the serial GA).
_ABSORB_SKIP_IND = {"peakArr", "bkgnArr", "new_bkgns", "shirley",
                    "solo_shirley", "peakDict", "pars_range"}


def _absorb_state(dst_ind, src_ind):
    """Copy a worker-evaluated individual's state back onto the original
    main-process object, preserving object identity (and therefore every
    aliasing behavior of the serial path). Bound-method attributes
    (peak.func, background.bkgn) and shared structures are kept."""
    for dst_p, src_p in zip(dst_ind.peakArr, src_ind.peakArr):
        for k, v in src_p.__dict__.items():
            if k not in ("func", "paramRange"):
                dst_p.__dict__[k] = v
    for name in ("bkgnArr", "new_bkgns", "shirley", "solo_shirley"):
        dst_list = getattr(dst_ind, name, None)
        src_list = getattr(src_ind, name, None)
        if dst_list is None or src_list is None:
            continue
        for dst_b, src_b in zip(dst_list, src_list):
            for k, v in src_b.__dict__.items():
                if k not in ("bkgn", "paramRange"):
                    dst_b.__dict__[k] = v
    for k, v in src_ind.__dict__.items():
        if k not in _ABSORB_SKIP_IND:
            dst_ind.__dict__[k] = v

if __name__ == "__main__":
    main()
