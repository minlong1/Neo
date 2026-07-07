import sys,glob,re,os
from tokenize import String
import numpy as np
import fnmatch
from PhysicsModules.XPS.xps_neo import xps_data
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.figure import Figure
import copy
from PhysicsModules.XPS.xps_neo.xps_individual import Individual
from PhysicsModules.XPS.xps_neo.xps_fit import peak, background
from PhysicsModules.XPS.xps_neo.loss import compute_loss
import scipy as scipy

"""
Author: Evan Restuccia evan@restuccias.com
"""


def sort_fold_list(dirs):
    fold_list = list_dirs(dirs)
    fold_list.sort(key=natural_keys)
    return fold_list
## Human Sort
def list_dirs(path):
    return [os.path.basename(x) for x in filter(
        os.path.isdir, glob.glob(os.path.join(path, '*')))]


def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]




class xps_analysis():
    """
    Contains all functions for the reading of and analysis of output files from the actualy genetic algorithm
    """
    def __init__(self,dirs,params, peakType, data_KE = False, data_XES = False, scale_var = False):
        """
        Params(Dictionary) : an dictionary of parameters

        Note, self.bkgns and self.peaks and lists of options for the possible types of peaks vs. backgrounds. This will allow the
        file reading portion to tell what is a background and what is a peak

        self.data_obj is an instance of the XPS_Data class
        self.is_singlet will become an array of bools where true = singlet, false = doublet

        """
        self.dirs = dirs
        fileName = params['fileName']
        self.bkgns = ['SVSC', 'Shirley', 'Baseline','Exponential', 'Linear', 'Polynomial 1', 'Polynomial 2', 'Polynomial 3', '3-Param Tougaard', '2-Param Tougaard']
        self.peaks = params['peaks']
        self.peak_options = self.peaks
        self.peakType = peakType
       
        for i in range(len(self.peaks)):
            self.peak_options[i] = self.peak_options[i].lower()
        self.bkgn_options = self.bkgns
        #self.numPeaks = len(self.peaks)
        for i in range(len(self.bkgns)):
            self.bkgns[i] = self.bkgns[i].lower()

        for i in range(len(self.peaks)):
            self.peaks[i] = self.peaks[i].lower()


        self.data_obj = params['data obj']
        self.is_singlet = []
        self.is_coster_kronig = []


        self.x = self.data_obj.get_x(data_KE, data_XES)
        self.y = self.data_obj.get_y(scale_var)
       

        #number of parameters that the GA fit
        #self.number_of_parameters = self.calculate_number_of_params(self.peaks,self.bkgns)

    def get_param_num(self,type):
        """
        (helper) takes in a peak/background type and matches it with the number of parameters it has
        always has the number of parameters + 1, the +1 is for it's type

        """
        
        #Peak Types:
        if type == 'voigt':
            return 5
        if type == 'gaussian':
            return 4
        if type == 'lorentzian':
            return 4
        if type == 'double lorentzian':
            return 6
        if type == 'doniach-sunjic':
            return 6

        #Backgrounds:
        if type == 'svsc':
            return 2
        if type == 'shirley':
            return 2
        if type == 'exponential':
            return 1
        if type == 'linear':
            return 3
        if type == 'baseline':
            return 2
        if type == 'polynomial 1':
            return 1
        if type == 'polynomial 2':
            return 1
        if type == 'polynomial 3':
            return 1
        if type == '3-param tougaard':
            return 1
        if type == '2-param tougaard':
            return 1
        else:
            print("Warning, need to add parameter number to get_param_num in analysis file")
            print("Type found was:" + type)
        
    def calculate_number_of_params(self,types):
        """
        peaks(array) array containing strings which represent the type of each peak, length is equal to the number of peaks
        bkgns(array) length is number of backgrounds, each element is a string representing the type of background

        (helper)
        Adds the number of parameters each type has together so the analysis will know how many columns to have in arrays of the parameters
        """

        num = 0

        for i,one_type in enumerate(types):

           
            num += self.get_param_num(one_type.lower().strip())

        return num

    def extract_data(self,scale_var,plot_err=True,passin=''):
        
        scale_var = scale_var
        
        """
        Extract data value using array data

        plot_err(bool) whether or not to plot the error
        passin(num?) Optional- the best fit

        After extracting data, evaluates all individuals and assigns self.best_ind to the best chisqr individual
        """
        #gets the parameters and best fit if they are not passed in
        if passin == '':
            full_param_list = []
            files = []
            for r, d, f in os.walk(self.dirs):
                f.sort(key = natural_keys)
                for file in f:
                    if fnmatch.fnmatch(file,'*_data.csv'):
                        files.append(os.path.join(r, file))
            files.sort(key=natural_keys)

            for i in range(len(files)):
                final_values,self.types = self.extract_from_result_file(files[i])
                full_param_list.append(final_values)
                if i == 0:
                    stack = full_param_list[i]
                stack = np.vstack((stack,full_param_list[i]))


            try: 
                self.types[0]
            except AttributeError:
                print("Analysis folder not selected")
                pass

            for i in range(len(self.types)): #Try selecting an output folder
                self.types[i]= self.types[i].strip()
            self.number_of_parameters = self.calculate_number_of_params(self.types)
            #print(full_param_list)
            self.full_matrix = stack
            #print(stack)
            bestFit,err = self.construct_bestfit_err_mat(stack,scale_var,plot_err)

            best_Fit = np.reshape(copy.deepcopy(bestFit),(1, -1))

            #print(best_Fit)
        else:
            err = np.zeros([self.number_of_parameters,1])
            err = passin
            bestFit = np.zeros([self.number_of_parameters,1])
            best_Fit = passin

        self.bestFit = best_Fit.flatten()
        self.err = err.flatten()
        #should print error here later, but need to customize it to number of params
        self.err = list(self.err)
        new_err = []
        new_best = []
        for i in range(len(self.err)):
            new_err.append(round(self.err[i],2))
            new_best.append(round(self.bestFit[i],2))
        #print("Parameter Error values are: ", new_err)



        error_array = []
        peaks = []
        bkgns = []
        bkgn_len_mod = 0
        bkgn_len_mod = 0
        for i,Type in enumerate(self.types):
            if Type.lower().strip() in self.peaks:
                peaks.append(Type)
            if Type.lower().strip() in self.bkgns:
                bkgns.append(Type)

                if Type == 'Linear':
                    bkgn_len_mod +=1

                if Type == '3-Param Tougaard':
                    bkgn_len_mod +=1
        error_length = self.number_of_parameters - len(peaks) - len(bkgns)

        #print(new_err)

        #Since error values come in one large array, I created a series of if and for statements that split up the array into the seperate peaks and sums them into an array of arrays
        #The background errors go to their own array
        l = 0
       
        k = 0

        peakNum = 0
        for peak in peaks: #iterates for as many peaks there are
            #Fixed error in logic --> Should work for all peaks now 11/06/2023
            #Issues when using two different peak types: Need to fix this, if two different peaks change range of second peak to (l, k+1)
            error = []
            #print("PEAK", peak)

            #Issue with error values and coster-kronig: All go to zero? Need to fix
            if(peak.lower() == 'voigt'):
               
                k += 4




                if self.is_singlet[peakNum] == False:
                    k += 2
                for j in np.arange(l, k):
                    error.append(new_err[j])
              
                l = k
                peakNum += 1
            elif(peak.lower() == 'gaussian'):
               
                
                k += 3
               
                #print(k, l)


                if self.is_singlet[peakNum] == False:
                    k += 2 #added -2 in for issue --> Not sure why it gets too big. Logic error
                    #print("doublet", k, l)
               
                for j in np.arange(l, k):
                    #print(k, j)
                    #print("ERROR", new_err[j])
                    #print(k, j)
                    #print("ERROR", new_err[j])
                    error.append(new_err[j])
               
                l = k
                peakNum += 1

            elif(peak.lower() == 'lorentzian'):
                k += 3




                if self.is_singlet[peakNum] == False:
                    k += 2
               
                for j in np.arange(l, k):
                    error.append(new_err[j])

                l = k
                peakNum += 1

            elif(peak.lower() == 'double lorentzian'):
                
                k += 5
               
                


                if self.is_singlet[peakNum] == False:
                    k += 2
              
                for j in np.arange(l, k):
                    #print(k, j)

                    #print("ERROR", new_err[j])
                   
                    error.append(new_err[j])
                    if new_err[j] == 0.0:
                        pass

                l = k
                peakNum += 1

            elif(peak.lower() == 'doniach-sunjic'):
                k += 5


                if self.is_singlet[peakNum] == False:
                    k += 2
              
                for j in np.arange(l, k):
                    error.append(new_err[j])

                l = k
                peakNum += 1

            error_array.append(error)
            

        error_bkgns = []
        bkgn_len = len(new_err) - len(bkgns) - bkgn_len_mod
        bkgn_len = len(new_err) - len(bkgns) - bkgn_len_mod
        #print(new_err)
        for i in reversed(range(bkgn_len, len(new_err))):
                error_bkgns.append(new_err[i])
        error_bkgns = error_bkgns[::-1] #flip it to appear in the same order as the printed out peak attributes
        self.error_bkgns = error_bkgns
        self.error_array = error_array





        self.bestFit_mat = best_Fit

        # print("Not Normalized:")
        self.best_Fit_n = copy.deepcopy(best_Fit)
        #print(self.best_Fit_n)
        # self.best_Fit_n[:,0] -= self.scaler.min_
        # self.best_Fit_n[:,0] /= self.scaler.scale_
        # print(self.best_Fit_n)
        # self.x_model = np.arange
        self.best_ind = self.create_individual_from_bestfit(self.best_Fit_n[0],scale_var)
        print("peak attributes read are:")
        #print(self.best_Fit_n)
        m, b, BTou2, BTou3, parameters = self.best_ind.get_analysis_params(self.x,self.y,self.is_coster_kronig)
        print(parameters)
        for bkgn in bkgns:
            if(bkgn.lower() == '3-param tougaard'):
                print("3-Parameter Tougaard values: ","B = ", round(BTou3, 2), "C = 1000", "D = 13300")
            elif(bkgn.lower() == '2-param tougaard'):
                print("2-Parameter Tougaard values: ", "C = 1643")
            elif(bkgn.lower() == 'exponential'):
                print("Exponential values: ", "m = ", m, "b = ", b) #How to get these values from xps_fit??
        self.y_model,self.peak_components,self.bkgn_components = self.best_ind.getFitWithComponents(self.x,self.y, self.peakType, self.bkgns, self.is_coster_kronig)
        self.residual = [0]*len(self.x)
        for j in range(len(self.x)):
            self.residual[j] = ((self.y_model[j]- self.y[j]))

        #Moved print statement down so that errors appear after the attributes:
        print("Peak Error Values are: ", error_array)
        print("Background Error Values are: ", error_bkgns)
        '''
        self.y_model_components = self.best_ind.get_peaks()
        self.bkgn_components = self.best_ind.get_backgrounds()
        self.totalbackground = [0] * len(self.x)
        for i in range(len(self.bkgn_components)):
            self.bkgn_components[i] = self.bkgn_components[i].getY(self.x,self.y)
            for k in range(len(self.x)):
                self.totalbackground[k] += self.bkgn_components[i][k]

        for i in range(len(self.y_model_components)):
            self.y_model_components[i] = self.y_model_components[i].getY(self.x) + self.totalbackground
        '''

    def extract_from_result_file(self,file):
        """
        Pulls the data from each given file in the selected folder

        When reading these files, it first generates an array of all items in the CSV
        Then it will check the type of the object read.

        When finished, this function will have updated self.types, self.is_singlet, and return csv_numbers
        self.types = array where each element is the type of the next item in the file, i.e.
            1,2,3,4,true,0,'voigt',.5,'Shirley' would read as 'voigt','Shirley
        self.is_singlet is the same as self.types except it only has the booleans
        csv_numbers is all the numbers left over

        The reasoning for 3 different arrays is the original analysis uses reshape often, which
        requires all the same type

        """
        try:
            os.path.exists(file)
            csv_unflatten = []
            csv_numbers = []
            gen_csv = np.genfromtxt(file,delimiter=',',dtype = None,encoding= None)
            for row in gen_csv:
                temp_row = []
                temp_numbers = []
                types = []
                is_doublet = False
                peaks = []
                for k in range(len(row)):
                    typeCheck = type(row[k])

                    if k == len(row)-1 or k == len(row):
                        pass
                    else:
                        type_next = type(row[k+1])

                    if(typeCheck == np.float64):
                        temp_numbers.append(np.float64(row[k]))
                    elif(typeCheck == np.str_):
                        #types.append(row[k])
                        types.append(row[k])
                    elif(typeCheck == np.bool_) and (type_next == np.float64): #DOUBLET BOOL                        
                        self.is_singlet.append(row[k]) 
                        
                    elif(typeCheck == np.bool_) and (type_next == np.str_): #SINGLET OR IS_COSTER_KRONIG
                       type_before = type(row[k-3])
                       type_be4 = type(row[k-4])
                       type_b_4 = type(row[k-5])
                       

                       #Coster-Kronig Bool comes after spin-orbit splitting value..
                       if(type_before == np.bool_) or (type_be4 == np.bool_) or (type_b_4 == np.bool_): #Should be the new peak width value if coster-kronig
                           self.is_coster_kronig.append(row[k])
                       else: #Should be the amplitude if it is a singlet 
                            self.is_singlet.append(row[k]) 
                        
                    temp_row.append(row[k])


                if len(csv_unflatten) == 0:
                        csv_unflatten = temp_row
                        csv_numbers = temp_numbers
                else:
                    csv_unflatten = np.vstack((csv_unflatten,temp_row))
                    csv_numbers = np.vstack((csv_numbers,temp_numbers))

            #csv_unflatten = gen_csv.reshape((len(gen_csv),self.number_of_parameters+1))
            #print("Type is " + str(gen_csv[0]))
            return csv_numbers,types
        except OSError:
            print(" " + str(file) + " Missing")
            pass

    def  construct_bestfit_err_mat(self,full_matrix,scale_var,plot=False):
        """
        Constructs a matrix of all fits, where each row is an individual and the column is the parameters
        Gets the score of all of these indiiduals, and assigns the bestFit object from that, note that the bestFit is
        just the array of best parameters, not a full individual yet

        full_matrix <matrix>: n by m where n is the number of samples and m is the number of parameters
        npaths <int>: number of paths
        plot <bol>: if it plot or not

        """
        score = []
        full_mat_var_cov = np.cov(full_matrix.T)
        full_mat_diag = np.diag(full_mat_var_cov)
        err = np.sqrt(full_mat_diag)

        labels = self.generate_labels()
#         print(full_mat)
#         bestFit = np.mean(full_mat,axis=0)
#         bestFit = full_mat[-1,:]

        for i in self.full_matrix:
            #arr = i.reshape((1,-1))
            ind = self.create_individual_from_bestfit(i,scale_var)
            score.append(self.fitness(ind))
        bestScore = np.nanargmin(score)
        bestFit = full_matrix[bestScore]
        if plot:
            plt.figure(figsize=(6.4,4.8))
            plt.xticks(np.arange(len(full_mat_diag)),labels[0],rotation=70);
            plt.bar(np.arange(len(full_mat_diag)),np.sqrt(full_mat_diag))

        return bestFit,err

    def create_individual_from_bestfit(self, fit, scale_var):
        """
        fit(list) a list of parameters in order
        First the function creates an empty individual object
        The array uses the self.types list and the get_param_num to get the correct elements of the list,
        and then uses the set function to hand all of those values to each individual peak and/or background as well as with
        its doublet information if is_singlet indicates it is a doublet

        Returns the individual made up of the parameters from the fit
        """
        pars_range = '' #dummy range, since we're going to set the params again later

        #construct peaks and backgrounds list for individual constructor
        peaks = []
        bkgns = []
        j = 0
        count = 0 #Added in to help with multiple shirley backgrounds
        for i,Type in enumerate(self.types):
            if Type.lower().strip() in self.peaks:
                peaks.append(Type)
            if Type.lower().strip() in self.bkgns:
                bkgns.append(Type)
                
        
        #creates dummy individual, with random attributes
        #Phase 3d: the package Individual expects the deduplicated
        #background type list (per-peak SVSC objects are created inside);
        #the per-row list from the CSV drives the setter loop below.
        bkgns_unique = []
        for b in bkgns:
            if b not in bkgns_unique:
                bkgns_unique.append(b)
        ind = Individual(bkgns_unique,peaks,scale_var,pars_range)
        peaks2 = ind.get_peaks()
        bkgns2 = ind.get_backgrounds()
        
        #going to custom set the attributes of the individual to match the ones in the fit
        current_index = 0
        peak_index = 0
        bkgn_index = 0
        ck_index = 0
        


        for i,Type in enumerate(self.types):
            num_params = self.get_param_num(Type.lower())
            

            param_list = []
            doublet_params = []
            CK_params = []
            CK_index_2 = 1
            for k in range(current_index,current_index+num_params-1):
                param_list.append(fit[k]) #error here if singlet, then doublet, then singlet in peak order. Out of bounds error
                current_index = k+1

            if Type.lower() in self.peak_options:
                if not self.is_singlet[peak_index]:
                    for k in range(current_index, current_index+2):
                        doublet_params.append(fit[k])
                        current_index = k+1
                   
                    if self.is_coster_kronig[ck_index] == True:
                        if Type.lower() == 'double lorentzian':
                            CK_index_2 = 2
                        for j in range(current_index, current_index+CK_index_2):
                            CK_params.append(fit[j])
                            current_index = j + 1
                            
                        ck_index += 1
                    else:
                        ck_index += 1
                    
                if self.is_singlet[peak_index] == True:
                    temp_ck = False
                elif self.is_singlet[peak_index] == False:
                    temp_ck = self.is_coster_kronig[ck_index -1]
                else: #Backgrounds
                    temp_ck = False
                peaks2[peak_index].set(param_list,doublet_params, CK_params,self.is_singlet[i],temp_ck)
                peak_index += 1
            
         
            #At later point in time, we could move this to a single set function in xps fit like the peaks
            if Type.lower() == 'baseline':
                bkgns2[bkgn_index].set_baseline(param_list)
                bkgn_index+=1
                count += 1
            elif Type.lower() == 'svsc':
                
                value =  bkgns2[count].set_shirley_sherwood(param_list) #IDK why but this needs to be here otherwise shirley = -1 when coster-kronig selected 
                peak_num = len(peaks) - 1
                count+=1
            elif Type.lower() == 'shirley':

                value =  bkgns2[bkgn_index].set_shirley(param_list) #IDK why but this needs to be here otherwise shirley = -1 when coster-kronig selected 
                bkgn_index+=1
                count+=1
                '''
                if j <= peak_num: #This isnt being seen for some reason...
                    bkgns2[count].set_shirley_sherwood(param_list) #For some reason this line will throw an error when selecting both shirley and linear background type. Shirley has to be selected first for this error to go away      
                    j += 1
                    bkgn_index +=1 #Changed from 2 to 1 because I got rid of shirley background parameter since it wasnt doing anything.
                    count += 1

 
                else:
                    pass
                '''
                '''
                if j < len(peaks):
                    j += 1
                    pass
                else:
                    bkgn_index +=2 #changed from 1 to 2 for new Shirley Background using areas
                '''



            elif Type.lower() == 'linear':

                bkgns2[bkgn_index].set_linear(param_list)
                bkgn_index+=1
                count += 1

            elif Type.lower() == '3-param tougaard':

                bkgn_index+=1
                count += 1

            elif Type.lower() == '2-param tougaard':

                bkgn_index+=1 #Should this be commented out since it has no parameters in best fit array?
                count += 1


            elif Type.lower() == 'exponential':
                bkgn_index+=1
                count += 1

            elif Type.lower() == 'polynomial 1':

                bkgn_index+=1
                count += 1

            elif Type.lower() == 'polynomial 2':

                bkgn_index+=1
                count += 1

            elif Type.lower() == 'polynomial 3':

                bkgn_index+=1
                count += 1

        return ind



    def generate_labels(self):
        """
        I don't know what this one does yet, I never used it
        """
        label=[]
        amp_label = []
        center_label = []
        sigma_label = []
        gamma_label = []
        asymmetry_label = []
        asymmetryDoniach_label = []

        for i in range(1,self.number_of_parameters+1):
            label.append('amp_' + str(i))
            amp_label.append('amp_' + str(i))

            label.append('center_' + str(i))
            center_label.append('center_' + str(i))

            label.append('sigma_' + str(i))
            sigma_label.append('sigma_' + str(i))

            label.append('gamma_' + str(i))
            gamma_label.append('gamma_' + str(i))

            label.append('asymmetry_' + str(i))
            asymmetry_label.append('asymmetry_' + str(i))

            label.append('asymmetryDoniach_' + str(i))
            asymmetryDoniach_label.append('asymmetryDoniach_' + str(i))

        return label,amp_label,center_label,sigma_label,gamma_label,asymmetry_label,asymmetryDoniach_label



    def FWHM(self, indObj):
        FWHM = [0]*len(self.peaks)
        Individual = indObj
        FWHM = Individual.getFWHM(self.x, self.peakType)

        return FWHM



    def fitness(self,indObj):
        """
        Evaluate fitness of a individual

        Phase 3 unification: the loss math is xps_neo.loss.compute_loss,
        the same function the GA optimizes. apply_penalties=False reports
        the GA's final-phase definition (genNum >= gen_alt), so the value
        shown here matches the "Best fit" the GA logs at the end of a run.
        Decision recorded in docs/fork-divergence.md; previous behavior
        (N=30 window, penalties always on) is in the analysis goldens'
        CHANGELOG entry.
        """
        Individual = indObj
        N = 30
        y_left = self.y[:N]
        y_right = self.y[-N:]
        self.y_left_avg = sum(y_left)/N
        self.y_right_avg = sum(y_right)/N

        yTotal = Individual.getFit(self.x, self.y, Individual.backgrounds)

        loss, _residual = compute_loss(self.x, self.y, yTotal,
                                       self.y_left_avg, self.y_right_avg,
                                       apply_penalties=False)
        return loss

    def score(self,verbose=False):
        
        loss =self.fitness(self.best_ind)
        print('Fitness Score (Chi2):', loss)
        # print('Fitness Score (ChiR2):', loss/(len(self.x_raw)-4*self.npaths))
        # self.fwhm(verbose=verbose)
        # self.cal_area(verbose=verbose)
        # print(self.err_full)

    def analyze(self):
        
        #should calculate area and other params
        #print("Full analysis not yet functional")

        def area(yVals):
            """
            (helper) calculate area as Trapezoidal Riemann sum
            """
            area = 0
            #y_total = len(yVals)*[0]
            for i in range(len(yVals)-1):
                """
                (yVals[i]  -  Background[i])  +  (yVals[i + 1]  +  Background[i + 1])
                ----------------------------------------------------------------------
                                (2 * (x[i + 1] - x[i])

                """
                area += (((yVals[i]-self.background[i])+(yVals[i+1]-self.background[i+1]))/2)*(self.x[i+1]-self.x[i])


            right_area = 0
            left_area = 0
           
            h = (self.x[-1] -self.x[0])/(len(yVals)-1)
            for i in range(1, len(yVals)):
                right_area += (yVals[i]-self.background[i])*h #Finding area with right endpoint riemann sum (positive/negative area error)
            for i in range(1, len(yVals)):
                left_area += (yVals[i-1]-self.background[i-1])*h #Finding area with left endpoint riemann sum (negative/positive area error)
          
            if left_area < area:
                lower_error = left_area - area
                upper_error = right_area - area
            else:
                lower_error = right_area - area
                upper_error = left_area - area

                #y_total[i] = yVals[i]-self.background[i]
            #new_area = scipy.integrate.simps(y_total, self.x) #Calculating area using simpson's rule
            #print("OLD", area, "NEW", new_area, "DIFF", area-new_area) # minimal difference between Simpson's rule and trapozoid rule
            return abs(area), lower_error, upper_error

        self.background = [0] * len(self.x)
        for i in range(len(self.bkgn_components)):
                for j in range(len(self.x)):
                    self.background[j] += self.bkgn_components[i][j]

        self.peak_component_areas = [0] * len(self.peak_components)
        self.peak_upper_errors = [0] * len(self.peak_components)
        self.peak_lower_errors = [0] * len(self.peak_components)
        k = 0
        for peak in self.peak_components:
            self.peak_component_areas[k], self.peak_lower_errors[k], self.peak_upper_errors[k] = area(peak)
            k += 1
       
        self.totalArea, self.lower_error, self.upper_error = area(self.y_model)
        '''
        print("Areas")
        print("Total Area Upper and Lower Error:", self.upper_error, self.lower_error)
        print("Total Area: ", self.totalArea)
        print("Individual Peak Area: ", self.peak_component_areas)
        print("Individual Peak Area Upper Errors: ", self.peak_upper_errors)
        print("Individual Peak Area Lower Errors: ", self.peak_lower_errors)
        '''
        return self.totalArea, self.peak_component_areas
        #What are the error values for the area? Are they correlated to the background error values???
        '''
        def FWHM_values(yVals):
            MAD = np.mean(np.absolute(yVals-np.mean(self.y))) #Mean Absolute Deviaition
            self.fwhm_g = 2*self.gaussian*np.sqrt(2*np.log(2))
            self.fwhm_l = 2*MAD
            self.fhwm_v = fwhm_l/2 + np.sqrt(pow(fwhm_l, 2)/4 + pow(fwhm_g, 2))
            return self.fwhm_v, self.fwhm_l, self.fwhm_g
        self.FWHM = FWHM_values(self.y_model)
        print(FWHM)
        #FWHM_values = peak.getFWHM(self)
        '''







    def plot_data(self,scale_var,data_KE, data_XES,title='Test',fig_gui=None):
        y_val = []
        if fig_gui == None:
        
            first = self.y[0]
            if first < 1:
                scale_val = round(1/first)*100
            else:
                scale_val = first

            if scale_var == True:
            
                y_val = self.y/first #Dividing every element by the first value 
            
                y_val = self.y*scale_val #Multiply by 1000 to scale it
               
            else:
                y_val = self.y
            self.y = y_val
                

            
            

            plt.rc('xtick', labelsize='12')
            plt.rc('ytick', labelsize='12')
            # plt.rc('font',size=30)
            plt.rc('figure',autolayout=True)
            plt.rc('axes',titlesize=12,labelsize=12)
            # plt.rc('figure',figsize=(7.2, 4.45))
            plt.rc('axes',linewidth=1)

            plt.rcParams["font.family"] = "Times New Roman"
            # fig,ax = plt.subplots(1,1,figsize=(6,4.5))
            # ax.plot(self.x_raw,self.y_raw,'ko-',linewidth=1,label='Data')
            plt.plot(self.x,self.y,'b--',linewidth=1,label='Data')
            # ax.plot(self.x_slice,self.y_slice,'r--',linewidth=1.2,label='Slice Data')
            plt.plot(self.x,self.y_model,'k',linewidth=1,label='Fit')
            #plt.plot(self.x_linear,self.y_linear,'--',color='tab:purple')
            plt.legend()
        else:
            #ax  = fig_gui.add_subplot(111)
            y_val = []
            first = self.y[0]
            if first < 1:
                scale_val = round(1/first)*100
            else:
                scale_val = first
        
           
            if scale_var == True:
                y_val = self.y/first #Dividing every element by the first value 
            
                y_val = self.y*scale_val #Multiply by 1000 to scale it
              
            else:
                y_val = self.y
            self.y = y_val
          
           #added in gridspec to help make residual plot
            gs = gridspec.GridSpec(2,1, height_ratios=[1,0.1])
            ax = fig_gui.add_subplot(gs[0])
            ax2 = fig_gui.add_subplot(gs[1])
            gs.update(hspace=0)
            ax.plot(self.x,self.y,color='gray', linestyle='None', marker='o', markeredgewidth='0.3', markerfacecolor='None', markeredgecolor='gray', markersize='5',label='Data')
            ax.plot(self.x,self.y_model,'k-',linewidth=1.2,label='Fit')
            ax2.plot(self.x,self.residual,'k-',linewidth=1.2)
            #plot peak componenets
            if(len(self.peak_components)>1):
                for i,peak in enumerate(self.peak_components):

                    ax.plot(self.x,peak,'-',linewidth = 1,label=('Peak ' + str(i+1)))
                    '''
                    if data_KE == False:
                        ax.plot(self.x,peak,'-',linewidth = 1,label=('Peak ' + str(i+1)))
                    else:
                        ax.plot(self.x[::-1],peak,'-',linewidth = 1,label=('Peak ' + str(i+1)))
                    '''

            self.background = [0] * len(self.x)
            for i in range(len(self.bkgn_components)):
                for j in range(len(self.x)):
                    self.background[j] += self.bkgn_components[i][j]
            ax.plot(self.x,self.background,color='dimgray',linestyle='dashed',linewidth =1, label='background')
           
            if data_KE == False:
                ax.invert_xaxis()
                ax2.invert_xaxis()
            else:
                pass
            if data_XES == True:
                ax.invert_xaxis()
                ax2.invert_xaxis()


            #ax.invert_xaxis()
            #ax2.invert_xaxis()
            #ax.plot(self.x_linear,self.y_linear,'--',color='tab:purple')
            ax.legend()
            #ax.set_xlabel("Binding Energy (eV)")
            ax.set_ylabel("Intensity (a.u.)")
            ax2.set_ylabel("Res.")
            if data_KE == False:
                ax2.set_xlabel("Binding Energy (eV)")
            else:
                ax2.set_xlabel("Kinetic Energy (eV)")

            if data_XES == True:
                ax2.set_xlabel("Photon Energy (eV)")
            #ax2.set_xlabel("Binding Energy (eV)")
            fig_gui.tight_layout()

    def get_unique(self, scale_var):

        self.best_ind = self.create_individual_from_bestfit(self.best_Fit_n[0],scale_var)

        self.best_ind.uniqueness(self.x,self.y, self.peakType, self.bkgns, self.is_coster_kronig)


    def get_params(self):
        #return self.params
        is_coster_kronig = self.is_coster_kronig
        m, b, BTou2, BTou3, parameters = self.best_ind.get_analysis_params(self.x,self.y,is_coster_kronig)


        


        errors = self.error_array
        errors_bkgns = self.error_bkgns
        upper_errors_peak_area = self.peak_upper_errors
        lower_errors_peak_area = self.peak_lower_errors
        peak_areas = self.peak_component_areas
        FWHM_values = self.FWHM(self.best_ind)
        peak_y_vals = self.peak_components
        totalFit = self.y_model
        background_fit = self.background
        
        residual = self.residual

        print("Areas")
        print("Total Area Upper and Lower Error:", self.upper_error, self.lower_error)
        print("Total Area: ", self.totalArea)
        print("Individual Peak Area: ", self.peak_component_areas)
        print("Individual Peak Area Upper Errors: ", self.peak_upper_errors)
        print("Individual Peak Area Lower Errors: ", self.peak_lower_errors)


        '''
        
        X_r = [N_ind / (N_ind - N_p)*N]*diff

        N_ind = number of data points in array 
        N_p number of fitting parameters
        
        '''
        y_raw = self.y
        diff = 0
        for i in range(len(self.x)):
            diff += pow((totalFit[i] - y_raw[i]),2) 
        res = 1.0
        #res = 0.48 #resolution of Kratis AXIS NOVA?
        #res = 0.48-1.42 #Phi 5600 resolution
        step_size = self.x[2] - self.x[1]
        res_step = res/step_size
        N = len(self.x) 
        N_ind = N/res_step #len(self.x)
        N_p = len(errors[0])*len(peak_areas) + len(errors_bkgns)
        chi2_r = (N_ind/((N_ind - N_p)*N))*diff
        print("Reduced Chi2: ", round(chi2_r,2))
        return parameters,errors, errors_bkgns, BTou2, BTou3, peak_areas, FWHM_values, peak_y_vals, totalFit, background_fit, residual, y_raw, upper_errors_peak_area, lower_errors_peak_area
