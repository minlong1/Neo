"""
Created 7/5/23
Created the file, left out set peak to manually set a peak and its
parameters, it seems non essential for now
-Evan Restuccia (evan@restuccias.com)
"""
from PhysicsModules.XPS.xps_neo.xps_fit import peak,background
import numpy as np
import sys

class Individual:
    def __init__(self,backgrounds,peaks,scale_var,pars_range=''):
        """
        backgrounds (array) Array where each element is the name of the background type desired
        peaks (array) array where each element is the name of the desired peakType
        pars_range (dict) each key is the name of a parameter with a tuple that contains the range the parameter is allowed to explore
        """

        self.scale_var = scale_var
        self.pars_range = pars_range

        if pars_range == '':
            # Analysis container mode (Phase 3d): parameters are set
            # explicitly from a fit result; no ranges, no RNG draws.
            self._init_analysis_container(backgrounds, peaks)
            return

        #both peaks and backgrounds are arrays of strings which represent the type of the background of the peak/bkgn
        self.nPeaks = len(peaks)

        self.nBackgrounds = len(backgrounds)
        self.backgrounds = backgrounds
        #self.peakArr = [None]* self.nPeaks
        self.bkgnArr = [None] * self.nBackgrounds
        self.peakDict = {}
        is_singlet = pars_range['is_singlet']
        is_coster_kronig = pars_range['is_coster_kronig']
        self.backgrounds = backgrounds
        """
        the Binding Energy needs to be personalized
        we take the range, which right now is something like (0,.2), i.e. the allowed variance in BE
        Then we add our binding energy to it to move the range to the right spot
        """
        amp_inputs = pars_range['Amplitude'] #Declaring here so that pars_range['Amplitude'] can become a 1d array in for loop. This method isnt the best but it works
        BE_inputs = pars_range['Binding Energy']
        asym_inputs = pars_range['Asymmetry']
        sos_inputs = pars_range['spinOrbitSplitting']
        br_inputs = pars_range['br_range']
    
        self.peakArr = [None]*len(peaks)
        sigma_inputs = pars_range['Gaussian']
        gamma_inputs = pars_range['Lorentzian']
        gamma_CK_inputs = pars_range['Lorentzian Coster-Kronig']
        asym_CK_inputs = pars_range['Asymmetry Range Coster-Kronig']

        #Need full and temp values for these parameters because their ranges are declared by the user to be +/- some value. Below we go through and add the actual value to this range
        #Example: BE = 285 BE_range = (-0.2, 0.2), change to make BE_range = (284.8, 285.2)
        self.full_amp_range = pars_range['Amplitude']
        self.full_sos_range = pars_range['spinOrbitSplitting']
        self.full_br_range = pars_range['br_range']
        self.full_asym_range = pars_range['Asymmetry']
        self.full_sigma_range = pars_range['Gaussian']
        self.full_gamma_range = pars_range['Lorentzian']
        self.full_gamma_CK_range = pars_range['Lorentzian Coster-Kronig']
        self.full_asym_CK_range = pars_range['Asymmetry Range Coster-Kronig']
        self.full_energy_range = pars_range['Binding Energy']
    
        #print("FULL RANGE GAMMA CK", self.full_gamma_CK_range)
        

        num_peaks = pars_range['npeaks']
        
       
        for i in range(self.nPeaks): 
           
            '''
            range_key = 'Binding Energy'
            guess_key = 'BE'
            BE_range = pars_range[range_key]
            BE1,BE2 =BE_range[0],BE_range[1]
            binding_energy = pars_range[guess_key][i]
            pars_range[range_key][0],pars_range[range_key][1] = [BE_range[0] + binding_energy, BE_range[1] + binding_energy,]
            '''

            #Want to add in correlation --> make BE_value/sigma_value/gamma_value/amp_value equal to the peak correlated or to multiplication of that value by something?
            #Have range become zero for that value?
            #Not sure how this will work with the random choice value used in xps_fit --> Will the correlation only follow the entry value?
            self.BE_correlated = int(pars_range['BE_correlated'][i])
           
            self.BE_correlated_mult = pars_range['BE_correlated_mult'][i]
            self.sigma_correlated = int(pars_range['sigma_correlated'][i])
            self.sigma_correlated_mult = pars_range['sigma_correlated_mult'][i]
            self.gamma_correlated = int(pars_range['gamma_correlated'][i])
            self.gamma_correlated_mult = pars_range['gamma_correlated_mult'][i]
            self.amp_correlated = int(pars_range['amp_correlated'][i])
            self.amp_correlated_mult = pars_range['amp_correlated_mult'][i]
            self.asym_correlated = int(pars_range['asym_correlated'][i])
            self.asym_correlated_mult = pars_range['asym_correlated_mult'][i]
            self.sos_correlated = int(pars_range['sos_correlated'][i])
            self.sos_correlated_mult = pars_range['sos_correlated_mult'][i]
            self.br_correlated = int(pars_range['br_correlated'][i])
            self.br_correlated_mult = pars_range['br_correlated_mult'][i]
           
         

            if self.BE_correlated > num_peaks:
                print("Energy Peak Correlation out of range in Parameter Ranges tab")
                sys.exit(1)
            if self.sigma_correlated > num_peaks:
                print("Sigma Peak Correlation out of range in Parameter Ranges tab")
                sys.exit(1)
            if self.gamma_correlated > num_peaks:
                print("Gamma Peak Correlation out of range in Parameter Ranges tab")
                sys.exit(1)
            if self.amp_correlated > num_peaks:
                print("Amplitude Peak Correlation out of range in Parameter Ranges tab")
                sys.exit(1)
            if self.asym_correlated > num_peaks:
                print("Asymmetry Peak Correlation out of range in Parameter Ranges tab")
                sys.exit(1)
            if self.asym_correlated > num_peaks:
                print("Spin-Orbit Splitting Peak Correlation out of range in Parameter Ranges tab")
                sys.exit(1)
            if self.br_correlated > num_peaks:
                print("Branching Ratio Peak Correlation out of range in Parameter Ranges tab")
                sys.exit(1)

            range_key_BE = 'Binding Energy'
            guess_key_BE = 'BE'
            BE_inputs = np.array(BE_inputs)

            if BE_inputs.ndim > 1: #If it is a 2D array or not
                BE_range = BE_inputs[self.BE_correlated -1] #2D array  #Takes in same range as correlated peak # -1 because i starts at zero
            else:
                BE_range = BE_inputs #1D array


            BE_value = pars_range[guess_key_BE][self.BE_correlated-1] #Takes in same input value as correlated peak
            BE_value = BE_value*self.BE_correlated_mult
            if isinstance(BE_range, float): #Need this here because pars_range gets called later and is a 1D array.
                BE_range = BE_inputs
            
            #Check to make sure we do not have negative values
            if BE_range[0] + BE_value < 0:
                pars_range[range_key_BE] = [0, BE_range[1] +BE_value, BE_range[2]]
                
            else:
                pars_range[range_key_BE] = [BE_range[0] +BE_value, BE_range[1] +BE_value, BE_range[2]]

            self.temp_energy_range = pars_range[range_key_BE]

 
           
            #print("Calculated range is " + str(pars_range[range_key][0]) + " " + str(pars_range[range_key][1]))
            #pars_range['branching_ratio'] = pars_range['branching_ratios'][i]
            '''
            range_key_br = 'branching_ratio_range'
            guess_key_br = 'branching_ratios'
            br = pars_range[guess_key_br][i]
            
            branchingRatio_range = pars_range[range_key_br]
       
          
            BR1, BR2 =branchingRatio_range[0],branchingRatio_range[1]

            pars_range[range_key_br][0],pars_range[range_key_br][1] = [branchingRatio_range[0]+br, branchingRatio_range[1]+br,]
         
            #pars_range['spinOrbitSplit'] = pars_range['spinOrbitSplit'][i]
            
            range_key_sos = 'spinOrbitSplitting'
            guess_key_sos = 'SOS'
            spinOrbitSplit_range = pars_range[range_key_sos]
        
            SOS1,SOS2 =spinOrbitSplit_range[0],spinOrbitSplit_range[1]
            sos = pars_range[guess_key_sos][i]
         
            pars_range[range_key_sos][0],pars_range[range_key_sos][1] = [spinOrbitSplit_range[0] + sos, spinOrbitSplit_range[1] + sos,]
            '''

        


            #New Sigma ranges for each peak
            range_key_sigma = 'Gaussian'
            guess_key_sigma = 'Sigma'

            sigma_inputs = np.array(sigma_inputs)
            if sigma_inputs.ndim > 1:
                sigma_range = sigma_inputs[self.sigma_correlated-1] #2D array
            else:
                sigma_range = sigma_inputs #1D array

           
            sigma_value = pars_range[guess_key_sigma][self.sigma_correlated-1]
            sigma_value = sigma_value*self.sigma_correlated_mult
            if isinstance(sigma_range, float): #Need this here because pars_range gets called later and is a 1D array.
                sigma_range = sigma_inputs

            #Check to make sure we do not have negative values
            if sigma_range[0] + sigma_value < 0:
                pars_range[range_key_sigma] = [0, sigma_range[1] +sigma_value, sigma_range[2]]
            else:
                pars_range[range_key_sigma] = [sigma_range[0] +sigma_value, sigma_range[1] +sigma_value, sigma_range[2]]

            self.temp_sigma_range = pars_range[range_key_sigma]
           

            #New Gamma ranges for each peak
            range_key_gamma = 'Lorentzian'
            guess_key_gamma = 'Gamma'

            gamma_inputs = np.array(gamma_inputs)
            if gamma_inputs.ndim > 1:
                gamma_range = gamma_inputs[self.gamma_correlated-1] #2D array
            else:
                gamma_range = gamma_inputs #1D array

           
            gamma_value = pars_range[guess_key_gamma][self.gamma_correlated-1]
            gamma_value = gamma_value*self.gamma_correlated_mult
            if isinstance(gamma_range, float): #Need this here because pars_range gets called later and is a 1D array.
                gamma_range = gamma_inputs

            #Check to make sure we do not have negative values
            if gamma_range[0] + gamma_value < 0:
                pars_range[range_key_gamma] = [0, gamma_range[1] +gamma_value, gamma_range[2]]
            else:
                pars_range[range_key_gamma] = [gamma_range[0] +gamma_value, gamma_range[1] +gamma_value, gamma_range[2]]

            self.temp_gamma_range = pars_range[range_key_gamma]

            
           



            '''
            #New user inputs of Sigma, Gamma, and Amplitude
            range_key_sigma = 'Gaussian'
            guess_key_sigma = 'Sigma'
            sigma_range = pars_range[range_key_sigma]
            S1,S2 =sigma_range[0],sigma_range[1]
            sigma_value = pars_range[guess_key_sigma][i]
            pars_range[range_key_sigma][0],pars_range[range_key_sigma][1] = [sigma_range[0] + sigma_value, sigma_range[1] + sigma_value,]

            range_key_gamma = 'Lorentzian'
            guess_key_gamma = 'Gamma'
            fwhm_range = pars_range[range_key_gamma]

            G1,G2 =fwhm_range[0],fwhm_range[1]
            gamma_value = pars_range[guess_key_gamma][i]
            pars_range[range_key_gamma][0],pars_range[range_key_gamma][1] = [fwhm_range[0] + gamma_value, fwhm_range[1] + gamma_value,]
            '''



            #Coster-Kronig Lorentz values:
            range_key_gamma_CK = 'Lorentzian Coster-Kronig'
            guess_key_gamma_CK = 'Gamma Coster-Kronig'
            #fwhm_range_CK = pars_range[range_key_gamma_CK]
            
            #CK1,CK2 =fwhm_range_CK[0],fwhm_range_CK[1]
            gamma_value_CK = pars_range[guess_key_gamma_CK][i]
            #print("GAMMA VALUE CK", gamma_value_CK)
            #pars_range[range_key_gamma_CK][0],pars_range[range_key_gamma_CK][1] = [fwhm_range_CK[0] + gamma_value_CK, fwhm_range_CK[1] + gamma_value_CK]

            #self.temp_gamma_CK_range = pars_range[range_key_gamma_CK]
            #print("GAMMA CK RANGE", self.temp_gamma_CK_range)



            gamma_CK_inputs = np.array(gamma_CK_inputs)
            if gamma_CK_inputs.ndim > 1:
                gamma_CK_range = gamma_CK_inputs[self.gamma_correlated-1] #2D array
            else:
                gamma_CK_range = gamma_CK_inputs #1D array

           
            gamma_CK_value = pars_range[guess_key_gamma_CK][self.gamma_correlated-1]
            gamma_CK_value = gamma_CK_value*self.gamma_correlated_mult
            if isinstance(gamma_CK_range, float): #Need this here because pars_range gets called later and is a 1D array.
                gamma_CK_range = gamma_CK_inputs

            #Check to make sure we do not have negative values
            if gamma_range[0] + gamma_value < 0:
                pars_range[range_key_gamma_CK] = [0, gamma_range[1] +gamma_value*2.0, gamma_CK_range[2]]
            else:
                pars_range[range_key_gamma_CK] = [1.1*gamma_value, gamma_range[1] +gamma_value*4.0, gamma_CK_range[2]] #Changed range so that gamma_ck has a minimum value of gamma*2 for that peak and maximum of 10*gamma 
                #pars_range[range_key_gamma_CK] = [0.5624998, 0.56250001, gamma_CK_range[2]]
            self.temp_gamma_CK_range = pars_range[range_key_gamma_CK]











            #New Amplitude ranges for each peak
            range_key_amp = 'Amplitude'
            guess_key_amp = 'Amp'

            amp_inputs = np.array(amp_inputs)
            if amp_inputs.ndim > 1:
                amp_range = amp_inputs[self.amp_correlated-1] #2D array
            else:
                amp_range = amp_inputs #1D array

            #A1,A2,A3 =amp_range[0],amp_range[1]
            amp_value = pars_range[guess_key_amp][self.amp_correlated-1]
            amp_value = amp_value*self.amp_correlated_mult
            #print("AMP VAL", amp_value, "PEAK", i)
            #pars_range[range_key_amp][i][0],pars_range[range_key_amp][i][1],pars_range[range_key_amp][i][2] = [amp_range[0] + amp_value, amp_range[1] + amp_value, amp_range[2]]

            if isinstance(amp_range, float): #Need this here because pars_range gets called later and is a 1D array.
                amp_range = amp_inputs

            #Check to make sure we do not have negative values
            if amp_range[0] + amp_value < 0:
                pars_range[range_key_amp] = [0, amp_range[1] +amp_value, amp_range[2]]
            else:
                pars_range[range_key_amp] = [amp_range[0] +amp_value, amp_range[1] +amp_value, amp_range[2]]

            self.temp_amp_range = pars_range[range_key_amp]



            #New Asymmetry ranges for each peak
            range_key_asym = 'Asymmetry'
            guess_key_asym = 'Asym'

            asym_inputs = np.array(asym_inputs)
            if asym_inputs.ndim > 1:
                asym_range = asym_inputs[self.asym_correlated-1] #2D array
          
            else:
                asym_range = asym_inputs #1D array
     
         
            asym_value = pars_range[guess_key_asym][self.asym_correlated-1]
            asym_value = asym_value*self.asym_correlated_mult
          
            if isinstance(asym_range, float): #Need this here because pars_range gets called later and is a 1D array.
                asym_range = asym_inputs

            #Check to make sure we do not have negative values
            if asym_range[0] + asym_value < 0:
                pars_range[range_key_asym] = [0, asym_range[1] +asym_value, asym_range[2]]
            else:
                pars_range[range_key_asym] = [asym_range[0] , asym_range[1] , asym_range[2]]

            self.temp_asym_range = pars_range[range_key_asym]





            #Coster-Kronig Asymmetry values:
            range_key_asym_CK = 'Asymmetry Range Coster-Kronig'
            guess_key_asym_CK = 'Asymmetry Coster-Kronig'

            asym_value_CK = pars_range[guess_key_asym_CK][i]



            asym_CK_inputs = np.array(asym_CK_inputs)
            if asym_CK_inputs.ndim > 1:
                asym_CK_range = asym_CK_inputs[self.asym_correlated-1] #2D array
            else:
                asym_CK_range = asym_CK_inputs #1D array

           
            asym_CK_value = pars_range[guess_key_asym_CK][self.asym_correlated-1]
            asym_CK_value = asym_CK_value*self.asym_correlated_mult
            if isinstance(asym_CK_range, float): #Need this here because pars_range gets called later and is a 1D array.
                asym_CK_range = asym_CK_inputs

            #Check to make sure we do not have negative values
            if asym_range[0] + asym_value < 0:
                pars_range[range_key_asym_CK] = [0, asym_range[1] +asym_value*2.0, asym_CK_range[2]]
            else:
                pars_range[range_key_asym_CK] = [1.1*asym_value, asym_range[1] +asym_value*4.0, asym_CK_range[2]] 

            self.temp_asym_CK_range = pars_range[range_key_asym_CK]













            #New Spin-Orbit Splitting ranges for each peak
            range_key_sos = 'spinOrbitSplitting'
            guess_key_sos = 'SOS'

            sos_inputs = np.array(sos_inputs)
            if sos_inputs.ndim > 1:
                sos_range = sos_inputs[self.sos_correlated-1] #2D array
          
            else:
                sos_range = sos_inputs #1D array
     
            SOS1,SOS2 =sos_range[0],sos_range[1]
            sos_value = pars_range[guess_key_sos][self.sos_correlated-1]
            sos_value = sos_value*self.sos_correlated_mult
       
          
            if isinstance(sos_range, float): #Need this here because pars_range gets called later and is a 1D array.
                sos_range = sos_inputs

            #Check to make sure we do not have negative values
            if sos_range[0] + sos_value < 0:
                pars_range[range_key_sos] = [0, sos_range[1] +sos_value, sos_range[2]]
            else:
                pars_range[range_key_sos] = [sos_range[0] +sos_value, sos_range[1] +sos_value, sos_range[2]]

            self.temp_sos_range = pars_range[range_key_sos]






            #New Branching Ratio ranges for each peak
            range_key_br = 'br_range'
            guess_key_br = 'BR'

            br_inputs = np.array(br_inputs)
            if br_inputs.ndim > 1:
                br_range = br_inputs[self.br_correlated-1] #2D array
          
            else:
                br_range = br_inputs #1D array
     
            BR1,BR2 =br_range[0],br_range[1]
            br_value = pars_range[guess_key_br][self.br_correlated-1]
            br_value = br_value*self.br_correlated_mult
            
            
          
            if isinstance(br_range, float): #Need this here because pars_range gets called later and is a 1D array.
                br_range = br_inputs

            #Check to make sure we do not have negative values
            if br_range[0] + br_value < 0:
                pars_range[range_key_br] = [0, br_range[1] +br_value, br_range[2]]
            else:
                pars_range[range_key_br] = [br_range[0] +br_value, br_range[1] +br_value, br_range[2]]

            self.temp_br_range = pars_range[range_key_br]
          


            


  
            self.peakArr[i] = peak(pars_range,peaks[i],is_singlet = is_singlet[i], is_coster_kronig = is_coster_kronig[i], BE_correlated = pars_range['BE_correlated'][i], BE_correlated_mult = pars_range['BE_correlated_mult'][i], sigma_correlated = pars_range['sigma_correlated'][i],sigma_correlated_mult = pars_range['sigma_correlated_mult'][i],gamma_correlated = pars_range['gamma_correlated'][i],gamma_correlated_mult = pars_range['gamma_correlated_mult'][i],amp_correlated = pars_range['amp_correlated'][i],amp_correlated_mult = pars_range['amp_correlated_mult'][i],asym_correlated = pars_range['asym_correlated'][i],asym_correlated_mult = pars_range['asym_correlated_mult'][i],sos_correlated = pars_range['sos_correlated'][i],sos_correlated_mult = pars_range['sos_correlated_mult'][i], br_correlated = pars_range['br_correlated'][i],br_correlated_mult = pars_range['br_correlated_mult'][i]) #ERROR HERE --> HOW TO MAKE XPS_FIT ONLY SEE FIRST ARRAY IN AMPLITUDE???
            
            pars_range[range_key_sos][0],pars_range[range_key_sos][1] = SOS1,SOS2
            pars_range['br_range'][0],pars_range['br_range'][1] = BR1,BR2

        '''
        except:
            #Special option if youre creating a custom individual(i.e. for analysis)
            if pars_range =='':
                pass
            else:
                print("Error modding guesses")
                sys.exit(1)
        '''

        #each index in the peaks/background array is the name of the peak/background type to be used
        k=0
        n=0
        '''
        for i in range(self.nBackgrounds):
            if backgrounds[i] == 'SVSC':
                self.bkgnArr = [None] * (self.nBackgrounds + self.nPeaks)
        '''

        self.new_bkgns = [None]*self.nBackgrounds #new bkgn for making extra shirley backgrounds
        if 'SVSC' in backgrounds:
            self.new_bkgns = [None]*(len(peaks) + self.nBackgrounds -1) #Size of number of backgrounds + number of peaks -1 for overcounting shirley bkgn
            self._new_bkgns_names = []
        else:
            self._new_bkgns_names = []
        self.is_shirley = False
        self.is_solo_shirley = False
        j = 0
        for i in range(self.nBackgrounds):
            
            if backgrounds[i] == 'SVSC':

                #while n >= self.nPeaks:
                #self.bkgnArr[k] = background(pars_range,backgrounds[i], peaks[n]) #Creating a unique shirley background for each peak...I think

                #k +=1
                #n +=1
                self.shirley = [None]*len(peaks)
                for n in range(len(peaks)):
                    self._new_bkgns_names.append(backgrounds[i])
                    self.shirley[n] = background(pars_range,backgrounds[i], peaks[n])
                    self.new_bkgns[j] = background(pars_range,backgrounds[i], peaks[n])
                    j+=1
                    
                  
                self.is_shirley = True
                self.bkgnArr[k] = background(pars_range,'Baseline', peaks[n]) #Should be zero...Making baseline right now 
                k+=1

            elif backgrounds[i] == 'Shirley':   #Shirley issue: This is hard coded to only take the Shirley background of the first peak.  We want the overall fit instead. How to get overall peak area calc. and bkgn subtration when calculating Shirley peak?
               
                self.solo_shirley = [None]*len(peaks)
                self.is_solo_shirley = True
                self._new_bkgns_names.append(backgrounds[i])
                self.solo_shirley[0] = background(pars_range,backgrounds[i], peaks[n])
                self.new_bkgns[j] = background(pars_range,backgrounds[i], peaks[n])

                self.bkgnArr[k] = background(pars_range,'Baseline', peaks[n]) 
                j+=1
                k+=1
            else:
                #n = 0 #Doesn't matter which peak we are looking at --> Only matters in Shirley case
                self.bkgnArr[k] = background(pars_range,backgrounds[i], peaks[n])
                self.new_bkgns[j] = background(pars_range,backgrounds[i], peaks[n])
                self._new_bkgns_names.append(backgrounds[i])

                k+=1
                j+=1

        for i in range(self.nPeaks): 
          
            #pars_range[range_key][0],pars_range[range_key][1] = BE1,BE2
            for i in range(len(pars_range[range_key_BE])): pars_range[range_key_BE][i] = BE_range[i]
            #pars_range[range_key_sos][0],pars_range[range_key_sos][1] = SOS1,SOS2

            #Added in sigma, gamma, and amplitude
            for i in range(len(pars_range[range_key_sigma])): pars_range[range_key_sigma][i] = sigma_range[i]
            for i in range(len(pars_range[range_key_gamma])): pars_range[range_key_gamma][i] = gamma_range[i]
            #pars_range[range_key_sigma][0],pars_range[range_key_sigma][1] = S1,S2
            #pars_range[range_key_gamma][0],pars_range[range_key_gamma][1] = G1,G2
            for i in range(len(pars_range[range_key_gamma_CK])): pars_range[range_key_gamma_CK][i] = gamma_CK_range[i]
            #pars_range[range_key_gamma_CK][0],pars_range[range_key_gamma_CK][1] = CK1,CK2
            for i in range(len(pars_range[range_key_amp])): pars_range[range_key_amp][i] = amp_range[i]
            for i in range(len(pars_range[range_key_asym])): pars_range[range_key_asym][i] = asym_range[i]
            #pars_range[range_key_amp][0],pars_range[range_key_amp][1], pars_range[range_key_amp][2] = A1,A2,A3

        self.pars_range['Gaussian'] = self.full_sigma_range
        self.pars_range['Lorentzian'] = self.full_gamma_range
        self.pars_range['Lorentzian Coster-Kronig'] = self.full_gamma_CK_range
        self.pars_range['Asymmetry Range Coster-Kronig'] = self.full_asym_CK_range
        self.pars_range['Amplitude'] = self.full_amp_range
        self.pars_range['spinOrbitSplitting'] = self.full_sos_range
        self.pars_range['br_range'] = self.full_br_range
        self.pars_range['Asymmetry'] = self.full_asym_range
        self.pars_range['Binding Energy'] = self.full_energy_range
      
        

        '''
        for n in range(self.nPeaks):
            for i in range(self.nBackgrounds): #May have to iterate through peaks too if we want background to be peak dependent
                self.bkgnArr[k] = background(pars_range,backgrounds[i], peaks[n])
                k+=1
        '''
         # Create dictionary of peaks and backgrounds
         #MAYBE THE ISSUE IS HERE???
        for i in range(self.nPeaks):
            
            self.peakDict[f'peak_{i}'] = self.peakArr[i]

        for i in range(self.nBackgrounds):
            self.peakDict[f'bkgn_{i}'] = self.bkgnArr[i]
    

    def add_peak(self,peakType):
        self.peakArr.append(peak(self.pars_range,peakType))
    def add_bkgn(self,bkgnType):
        self.bkgnArr.append(background(self.pars_range,bkgnType))
    
   
    #adds all backgrounds and peaks as one y value array
    def getFit(self,x,y, backgrounds):

        shirley = []*self.nPeaks
        peakVals = []*self.nPeaks

        yFit = [0]*len(x)
        self._analysis_bkgn_parts = [] #Phase 3d: per-background curves for the analysis; appends only, no math change
        FWHM = [0]*self.nPeaks
        shirley = []*self.nPeaks
        solo_shirley = []
        peakVals = []*self.nPeaks
        for i in range(self.nPeaks):

            BE_correlated = int(self.peakArr[i].BE_correlated)
           
            BE_correlated_mult = self.peakArr[i].BE_correlated_mult
            sigma_correlated = int(self.peakArr[i].sigma_correlated)
            sigma_correlated_mult = self.peakArr[i].sigma_correlated_mult
            gamma_correlated = int(self.peakArr[i].gamma_correlated)
            gamma_correlated_mult = self.peakArr[i].gamma_correlated_mult
            amp_correlated = int(self.peakArr[i].amp_correlated)
            amp_correlated_mult = self.peakArr[i].amp_correlated_mult
            asym_correlated = int(self.peakArr[i].asym_correlated)
            asym_correlated_mult = self.peakArr[i].asym_correlated_mult
            sos_correlated = int(self.peakArr[i].sos_correlated)
            sos_correlated_mult = self.peakArr[i].sos_correlated_mult
            br_correlated = int(self.peakArr[i].br_correlated)
            br_correlated_mult = self.peakArr[i].br_correlated_mult
            
            asymD = self.peakArr[i].asymmetryDoniach
            BE = self.peakArr[BE_correlated-1].bindingEnergy
            BE = BE*BE_correlated_mult
            BE = round(BE, 5)
            width = self.peakArr[gamma_correlated-1].lorentz
            width = width*gamma_correlated_mult
            sigma = self.peakArr[sigma_correlated-1].gaussian
            sigma = sigma*sigma_correlated_mult
            A = self.peakArr[amp_correlated-1].amp
            A = A*amp_correlated_mult
            asym = self.peakArr[asym_correlated-1].asymmetry
            asym = asym*asym_correlated_mult
            #print("PEAK BE", BE)






            #asym = self.peakArr[i].asymmetry
            singlet = self.peakArr[i].is_singlet
            coster_kronig = self.peakArr[i].is_coster_kronig
            if coster_kronig == True:
                if self.peakArr[i].peakType.lower() == "double lorentzian":
                    asym_CK = round(self.peakArr[i].asym_CK, 7)
                else:
                    asym_CK = round(self.peakArr[asym_correlated-1].asymmetry,7)
                width_CK = round(self.peakArr[i].lorentz_CK, 7)
            else:
                width_CK = round(self.peakArr[gamma_correlated-1].lorentz,7)
                asym_CK = round(self.peakArr[asym_correlated-1].asymmetry,7)
            #branch = self.peakArr[i].branching_ratio
            #split = self.peakArr[i].spinOrbitSplit
            split = self.peakArr[sos_correlated-1].spinOrbitSplit
            split = split*sos_correlated_mult
            branch = self.peakArr[br_correlated-1].branching_ratio
            branch = branch*br_correlated_mult

     
         
          
            #print("PEAK #:", i)

            #Phase 4 (order-preserving/bit-exact): evaluate the peak curve once
            #and reuse it; the second getY call recomputed the identical array
            #(shape functions are deterministic and rebind peak_y each call).
            _peak_curve = self.peakArr[i].getY(x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split)
            yFit += _peak_curve
            peakVals.append(_peak_curve)
            FWHM += self.peakArr[i].getFWHM(x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split)
            
            #Shirley Background depends on number of peaks. Added to this section and removed from nBackgrounds below (dont need to get it twice)
            if self.is_shirley == True:
       
                self.pars_range['Gaussian'] = self.temp_sigma_range
                self.pars_range['Lorentzian'] = self.temp_gamma_range
                self.pars_range['Lorentzian Coster-Kronig'] = self.temp_gamma_CK_range
                self.pars_range['Amplitude'] = self.temp_amp_range
                self.pars_range['Asymmetry'] = self.temp_asym_range
                self.pars_range['spinOrbitSplitting'] = self.temp_sos_range
                self.pars_range['br_range'] = self.temp_br_range
                self.pars_range['Binding Energy'] = self.temp_energy_range
                #Get yFit --> use in Shirley for total area
                tot_area = 0
                shirley.append(self.shirley[i].getY(x,y, BE , width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split,self.scale_var, tot_area))
            

        if self.is_shirley == True:
            new_shirley = [0]*len(x)
            for n in range(self.nPeaks):
                for i in range(len(x)): #Adding all shirley bkgns together 
                    new_shirley[i] += shirley[n][i]
            #for i in range(len(x)):
            #    new_shirley[i] /= self.nPeaks #Making them the average of all of them 
           

            
                #new_shirley = np.divide(new_shirley, self.nPeaks) #dividing by the number of peaks to get the summed shirley background
            #for i in range(len(x)):
            #    print(x[i], new_shirley[i])
            yFit += new_shirley
            self._analysis_bkgn_parts.append(new_shirley)
        
        self.pars_range['Gaussian'] = self.full_sigma_range
        self.pars_range['Lorentzian'] = self.full_gamma_range
        self.pars_range['Lorentzian Coster-Kronig'] = self.full_gamma_CK_range
        self.pars_range['Asymmetry Range Coster-Kronig'] = self.full_asym_CK_range
        self.pars_range['Amplitude'] = self.full_amp_range
        self.pars_range['spinOrbitSplitting'] = self.full_sos_range
        self.pars_range['br_range'] = self.full_br_range
        self.pars_range['Asymmetry'] = self.full_asym_range
        self.pars_range['Binding Energy'] = self.full_energy_range
        
        tot_peak_fit = yFit
        self.tot_peak_fit = tot_peak_fit #Phase 3d: exposed for the analysis
        
        for i in range(self.nPeaks):
            #Shirley Background depends on number of peaks. Added to this section and removed from nBackgrounds below (dont need to get it twice)
            '''
            if self.is_shirley == True:
       
                self.pars_range['Gaussian'] = self.temp_sigma_range
                self.pars_range['Lorentzian'] = self.temp_gamma_range
                self.pars_range['Lorentzian Coster-Kronig'] = self.temp_gamma_CK_range
                self.pars_range['Amplitude'] = self.temp_amp_range
                self.pars_range['Asymmetry'] = self.temp_asym_range
                self.pars_range['spinOrbitSplitting'] = self.temp_sos_range
                self.pars_range['br_range'] = self.temp_br_range
                self.pars_range['Binding Energy'] = self.temp_energy_range
                #Get yFit --> use in Shirley for total area
                shirley.append(self.shirley[i].getY(x,y, BE , width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split,self.scale_var, tot_peak_fit))
            '''
            if self.is_solo_shirley == True:

                self.pars_range['Gaussian'] = self.temp_sigma_range
                self.pars_range['Lorentzian'] = self.temp_gamma_range
                self.pars_range['Lorentzian Coster-Kronig'] = self.temp_gamma_CK_range
                self.pars_range['Amplitude'] = self.temp_amp_range
                self.pars_range['Asymmetry'] = self.temp_asym_range
                self.pars_range['spinOrbitSplitting'] = self.temp_sos_range
                self.pars_range['br_range'] = self.temp_br_range
                self.pars_range['Binding Energy'] = self.temp_energy_range
                new_solo_shirley = [0]*len(x)
                if i >= 1:
                    pass
                else:
                    solo_shirley.append(self.solo_shirley[i].getY(x,y, BE , width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split,self.scale_var, tot_peak_fit))
                    for i in range(len(x)): #Adding all shirley bkgns together 
                        new_solo_shirley[i] += solo_shirley[0][i]
                    yFit += new_solo_shirley
                    self._analysis_bkgn_parts.append(new_solo_shirley)
            #for i in range(len(x)):
            #    new_shirley[i] /= self.nPeaks #Making them the average of all of them 
           

            
                #new_shirley = np.divide(new_shirley, self.nPeaks) #dividing by the number of peaks to get the summed shirley background
            #yFit += new_shirley
        
        self.pars_range['Gaussian'] = self.full_sigma_range
        self.pars_range['Lorentzian'] = self.full_gamma_range
        self.pars_range['Lorentzian Coster-Kronig'] = self.full_gamma_CK_range
        self.pars_range['Asymmetry Range Coster-Kronig'] = self.full_asym_CK_range
        self.pars_range['Amplitude'] = self.full_amp_range
        self.pars_range['spinOrbitSplitting'] = self.full_sos_range
        self.pars_range['br_range'] = self.full_br_range
        self.pars_range['Asymmetry'] = self.full_asym_range
        self.pars_range['Binding Energy'] = self.full_energy_range




        n = 0
        bkgn = []*self.nPeaks
        i = 0
        for o in range(self.nBackgrounds):
            
            BE_correlated = int(self.peakArr[n].BE_correlated)
            sigma_correlated = int(self.peakArr[n].sigma_correlated)
            gamma_correlated = int(self.peakArr[n].gamma_correlated)
            amp_correlated = int(self.peakArr[n].amp_correlated)
            asym_correlated = int(self.peakArr[n].asym_correlated)
            


            asymD = self.peakArr[n].asymmetryDoniach
            BE = self.peakArr[BE_correlated-1].bindingEnergy
            
            width = self.peakArr[gamma_correlated-1].lorentz
            sigma = self.peakArr[sigma_correlated-1].gaussian
            A = self.peakArr[amp_correlated-1].amp
            asym = self.peakArr[asym_correlated-1].asymmetry
            #asym = self.peakArr[n].asymmetry
            singlet = self.peakArr[n].is_singlet
            coster_kronig = self.peakArr[n].is_coster_kronig
            #branch = self.peakArr[n].branching_ratio
            #split = self.peakArr[n].spinOrbitSplit
            split = self.peakArr[sos_correlated-1].spinOrbitSplit
            split = split*sos_correlated_mult
            branch = self.peakArr[br_correlated-1].branching_ratio
            branch = branch*br_correlated_mult
 
           

            #print("BE", BE)
            #Need this here for calling Shirley so that it sees the proper ranges of each parameter
            if backgrounds[o] == 'SVSC':
                #Do nothing in this for loop because it is done in the nPeaks for loop
               
                for h in range(self.nPeaks):
                    '''
                    self.pars_range['Gaussian'] = self.temp_sigma_range
                    self.pars_range['Lorentzian'] = self.temp_gamma_range
                    self.pars_range['Lorentzian Coster-Kronig'] = self.temp_gamma_CK_range
                    self.pars_range['Amplitude'] = self.temp_amp_range
                    self.pars_range['Asymmetry'] = self.temp_asym_range
                    self.pars_range['spinOrbitSplitting'] = self.temp_sos_range
                    self.pars_range['br_range'] = self.temp_br_range
                    self.pars_range['Binding Energy'] = self.temp_energy_range
                    #yFit += self.new_bkgns[i].getY(x,y, BE , width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split,self.scale_var)
                    '''
                    i += 1
          
                
                n += 1 #Added this to loop through each peak but I believe we are only seeing it once
            elif backgrounds[o] == 'Shirley':
                for h in range(self.nPeaks):
                    i += 1
                n += 1
            else:
                _bkgn_contrib = self.new_bkgns[i].getY(x,y, BE , width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split,self.scale_var, tot_peak_fit) #This is where backgrounds are actually calculated --> Is Shirley taking into account the multiple peaks??
                yFit += _bkgn_contrib
                self._analysis_bkgn_parts.append(_bkgn_contrib)
                i +=1
           
        
        self.pars_range['Gaussian'] = self.full_sigma_range
        self.pars_range['Lorentzian'] = self.full_gamma_range
        self.pars_range['Lorentzian Coster-Kronig'] = self.full_gamma_CK_range
        self.pars_range['Asymmetry Range Coster-Kronig'] = self.full_asym_CK_range
        self.pars_range['Amplitude'] = self.full_amp_range
        self.pars_range['spinOrbitSplitting'] = self.full_sos_range
        self.pars_range['br_range'] = self.full_br_range
        self.pars_range['Asymmetry'] = self.full_asym_range
        self.pars_range['Binding Energy'] = self.full_energy_range
        
        #Turn peak ranges back into 2D array to sort through all peaks otherwise it defaults to the last peaks parameter ranges

        self.peak_components = peakVals #Phase 3d: for the analysis
        self.bkgn_components = self._analysis_bkgn_parts

        return yFit
    
    '''
    def defineRanges(self):
        self.pars_range['Gaussian'] = self.full_sigma_range
        self.pars_range['Lorentzian'] = self.full_gamma_range
        self.pars_range['Amplitude'] = self.full_amp_range
        self.pars_range['Asymmetry'] = self.full_asym_range
        self.pars_range['Binding Energy'] = self.full_energy_range
    '''
      
        
    def get(self):
        """
        Get the whole set
        """
        return (self.peakArr + self.bkgnArr)


    def get_params(self, for_output_file = False):
        params = []
        #fetches all the params as independent lists
       
        for i in range(len(self.peakArr)):

            params.append(self.peakArr[i].get(for_output_file)) #WHAT IS FOR_OUTPUT_FILE???
            shirley_params = []
          
            shirley_val = 0
        l = 0
       
        for h in range(len(self.backgrounds)): #Based on number of backgrounds ---> For Shirley we need to do this as many times as there are peaks.
            #params.append(self.bkgnArr[i].get())
         
            #This is for one k value --> We need one for each peak so they can stay and change with generations.
           
            if self._new_bkgns_names[l] == 'SVSC':
                for n in range(self.nPeaks):
                    shirley_params += self.shirley[n].get()
                   
                    for i in range(len(shirley_params)):
                        if (i % 2) == 0:
                            shirley_val = shirley_params[i] 
                        
                    #shirley_val = shirley_val/self.nPeaks #Need to get each k value not just the average....
    
                    shirley_param = [shirley_val, shirley_params[-1]]
                  
                    params.append(shirley_param)
                    l += 1
                 

            else:
                params.append(self.new_bkgns[l].get())
                l +=1
               
            

        #puts it in one array
        for i in range(1,len(params)):
            for k in range(len(params[i])):
                params[0].append(params[i][k])

        #print("Params : " + str(params[0]))
        return params[0]

    def get_peak(self,i):
        return self.peakArr[i].get()

    def get_peaks(self):
        return self.peakArr

    def get_background(self,i):
        return self.bkgnArr[i]
    def get_backgrounds(self):
        return self.bkgnArr

    def mutate_(self,chance):
        for peak in self.peakArr:
            peak.mutate(chance)
        for bkgn in self.bkgnArr:
            bkgn.mutate(chance)

    #forces a given peak to have the given values, returns 0 on success, -1 on failure
    def setPeak(self,i,param_arr):
        #param array comes in with its last element indicating its type
        peakType = param_arr[len(param_arr)-1]
        #if param_array is voigt, it comes in form [BE,Gauss,Lorentz,Amplitude,'Voigt']
        if peakType.lower() == 'voigt': #Idk if this is important or now???
            self.peakArr[i].set_voigt(param_arr)
            return 0
        elif peakType.lower() == 'double lorentzian':
            self.peakArr[i].set_doubleLorentz(param_arr)
            return 0
        elif peakType.lower() == 'gaussian':
            self.peakArr[i].set_gauss(param_arr)
            return 0
        elif peakType.lower() == 'lorentzian':
            self.peakArr[i].set_lorentz(param_arr)
            return 0
        elif peakType.lower() == 'doniach-sunjic':
            self.peakArr[i].set_doniachSunjic(param_arr)
            return 0
        else:
            
            return -1

    def setBkgn(self,i,param_arr): #I dont think this function is doing anything
        bkgnType = param_arr[len(param_arr)-1]
        #if param_array is voigt, it comes in form [BE,Gauss,Lorentz,Amplitude,'Voigt']
        if bkgnType.lower() == 'svsc':
            self.bkgnArr[i].set_shirley_sherwood(param_arr)

    def verbose(self):
        """
        Print out the Populations
        """
        for i in range(self.npaths):
            self.Population[i].verbose()
    
    def __len__(self):
        """Returns the length of the independent parameters

        Returns:
            int: length of list of parameters
        """
        return len(self.get_params())
    
    def checkBound(self):
        """
        Check if the parameters are within the bounds
        """

        count = 0
        for i in self.peakArr:
            #peak(self,peakType).checkOutbound(self)

            
            peak_num = count
            #Check to make sure this is the correct range for each peak not just the last peak
            self.pars_range['Gaussian'] = self.temp_sigma_range
            self.pars_range['Lorentzian'] = self.temp_gamma_range
            self.pars_range['Lorentzian Coster-Kronig'] = self.temp_gamma_CK_range
            self.pars_range['Amplitude'] = self.temp_amp_range
            self.pars_range['Asymmetry'] = self.temp_asym_range
            self.pars_range['spinOrbitSplitting'] = self.temp_sos_range
            self.pars_range['br_range'] = self.temp_br_range
            self.pars_range['Binding Energy'] = self.temp_energy_range

            i.checkOutbound(peak_num)
            count += 1
        self.pars_range['Gaussian'] = self.full_sigma_range
        self.pars_range['Lorentzian'] = self.full_gamma_range
        self.pars_range['Lorentzian Coster-Kronig'] = self.full_gamma_CK_range
        self.pars_range['Asymmetry Range Coster-Kronig'] = self.full_asym_CK_range
        self.pars_range['Amplitude'] = self.full_amp_range
        self.pars_range['spinOrbitSplitting'] = self.full_sos_range
        self.pars_range['br_range'] = self.full_br_range
        self.pars_range['Asymmetry'] = self.full_asym_range
        self.pars_range['Binding Energy'] = self.full_energy_range
        

    def __getitem__(self,i):
        """Gets the ith peaks

        Args:
            i (int): index of the parameter

        Returns:
            float: value of the parameter
        """
        return self.peakArr[i]


    def getPeaksType(self):
        """Gets the type of the peaks

        Returns:
            list: list of the peak types
        """
        return [i.peakType for i in self.peakArr]

    def setParams(self,params,paramType=list):
        """Sets the parameters of the individual

        Args:
            params (list): list of parameters
        """
        for i in self.peakDict:
            self.peakDict[i].set(params[i],paramType=paramType)

    #----------------- Phase 3d: analysis container mode -----------------
    # These methods let the post-analysis (gui/xps_analysis2.py) rebuild an
    # Individual from a fit result and evaluate it with the SAME getFit the
    # GA's fitness used, so displayed curves match what was optimized.

    # Phase 4.2: slim pickling for --workers. The transient evaluation
    # products (component curves) are large and regenerated by every
    # getFit call; the GA loop never reads them, so they are dropped from
    # worker round-trips. The main-process originals keep their values
    # (_absorb_state skips absent keys).
    _PICKLE_SKIP = ("peak_components", "bkgn_components",
                    "_analysis_bkgn_parts", "tot_peak_fit", "peakDict")

    def __getstate__(self):
        return {k: v for k, v in self.__dict__.items()
                if k not in self._PICKLE_SKIP}

    def _init_analysis_container(self, backgrounds, peaks):
        """Container construction (pars_range=''): no ranges, no RNG.

        `backgrounds` is the deduplicated type-name list (e.g.
        ['Baseline', 'SVSC']); per-peak SVSC objects are created to match
        the layout getFit expects (self.shirley / self.new_bkgns).
        Parameters are assigned afterwards via the set_* methods.
        """
        self.nPeaks = len(peaks)
        self.peaks = list(peaks)
        self.backgrounds = list(backgrounds)
        self.bkgn_names = self.backgrounds
        self.nBackgrounds = len(self.backgrounds)
        self.peakDict = {}
        self.pars_range = {}  # getFit's range bookkeeping writes junk keys here harmlessly
        self.BTou2 = 0
        self.BTou3 = 0
        self.tot_peak_fit = []
        self._analysis_container = True
        for attr in ("full_amp_range", "full_sos_range", "full_br_range",
                     "full_asym_range", "full_sigma_range", "full_gamma_range",
                     "full_gamma_CK_range", "full_asym_CK_range",
                     "full_energy_range", "temp_amp_range", "temp_sos_range",
                     "temp_br_range", "temp_asym_range", "temp_sigma_range",
                     "temp_gamma_range", "temp_gamma_CK_range",
                     "temp_asym_CK_range", "temp_energy_range"):
            setattr(self, attr, '')
        self.peakArr = [
            peak('', peaks[i], is_singlet=True, is_coster_kronig=False,
                 BE_correlated=i + 1, BE_correlated_mult=1,
                 sigma_correlated=i + 1, sigma_correlated_mult=1,
                 gamma_correlated=i + 1, gamma_correlated_mult=1,
                 amp_correlated=i + 1, amp_correlated_mult=1,
                 asym_correlated=i + 1, asym_correlated_mult=1,
                 sos_correlated=i + 1, sos_correlated_mult=1,
                 br_correlated=i + 1, br_correlated_mult=1)
            for i in range(self.nPeaks)
        ]
        self.new_bkgns = []
        self.bkgnArr = []      # per-fit-row objects, for serialization
        self._bkgn_rows = []   # per-fit-row type names
        self.is_shirley = False
        self.is_solo_shirley = False
        for name in self.backgrounds:
            if name == 'SVSC':
                self.shirley = [background('', 'SVSC', peaks[n])
                                for n in range(self.nPeaks)]
                self.new_bkgns.extend(self.shirley)  # placeholders; getFit's SVSC branch skips them positionally
                self.bkgnArr.extend(self.shirley)
                self._bkgn_rows.extend(['SVSC'] * self.nPeaks)
                self.is_shirley = True
            elif name == 'Shirley':
                self.solo_shirley = [background('', 'Shirley', peaks[0])]
                self.new_bkgns.append(self.solo_shirley[0])
                self.bkgnArr.append(self.solo_shirley[0])
                self._bkgn_rows.append('Shirley')
                self.is_solo_shirley = True
            else:
                obj = background('', name, peaks[0])
                self.new_bkgns.append(obj)
                self.bkgnArr.append(obj)
                self._bkgn_rows.append(name)
        for i in range(self.nPeaks):
            self.peakDict[f'peak_{i}'] = self.peakArr[i]
        for i in range(len(self.bkgnArr)):
            self.peakDict[f'bkgn_{i}'] = self.bkgnArr[i]

    def getFitWithComponents(self, x, y, peakType=None, backgrounds=None,
                             is_coster_kronig=None):
        """y_model plus per-peak and per-background component curves,
        computed by the same getFit the GA's fitness evaluates. The extra
        arguments are accepted for gui-era call compatibility and ignored.
        """
        yFit = self.getFit(x, y, self.backgrounds)
        # The analysis expects each peak component to sit on top of the
        # total background (it subtracts the background before computing
        # areas), matching the gui fork's composition.
        peak_components = []
        for comp in self.peak_components:
            total = np.asarray(comp, dtype=float).copy()
            for b in self.bkgn_components:
                total = total + np.asarray(b, dtype=float)
            peak_components.append(total)
        return yFit, peak_components, self.bkgn_components

    def get_analysis_params(self, x, y, is_coster_kronig):
        """Serialize the fitted parameters for the post-analysis (port of
        the gui fork's get_params(x, y, is_coster_kronig); the gui-global
        exponential m/b extraction had no package equivalent and reports
        0, 0 - see docs/fork-divergence.md)."""
        params = []
        ck_val = []
        temp_ck = 0
        for i in range(self.nPeaks):
            singlet = self.peakArr[i].is_singlet
            if singlet == True:
                ck_val.append(self.peakArr[temp_ck].is_coster_kronig)
                temp_ck += 1
            else:
                ck_val.append(False)

        for i in range(len(self.peakArr)):
            p = self.peakArr[i]
            saved_ck = p.is_coster_kronig
            p.is_coster_kronig = ck_val[i]
            params.append(p.get())
            p.is_coster_kronig = saved_ck

        exp_m, exp_b = 0, 0
        BTou2, BTou3 = 0, 0
        for i in range(len(self.bkgnArr)):
            params.append(self.bkgnArr[i].get())
            row = self._bkgn_rows[i].lower()
            if row == '2-param tougaard':
                BTou2 = self.bkgnArr[i].get_2Param_vals(x, y, self.tot_peak_fit)
            elif row == '3-param tougaard':
                BTou3 = self.bkgnArr[i].get_3Param_vals(x, y, self.tot_peak_fit)
        return exp_m, exp_b, BTou2, BTou3, params

    def getFWHM(self,x, peakType): #not seeing coster-kronig but doesnt matter here. Should we calculate FWHM for coster-kronig peak? What about asymmetric?
            FWHM = [0]*self.nPeaks
            for i in range(self.nPeaks):
                peakType = self.peakArr[i].peakType
                if(peakType == "Doniach-Sunjic"):
                    asymD = self.peakArr[i].asymmetryDoniach
                else:
                    asymD = 0
                BE = self.peakArr[i].bindingEnergy
                width = self.peakArr[i].lorentz
                sigma = self.peakArr[i].gaussian
                A = self.peakArr[i].amp
                asym = self.peakArr[i].asymmetry
                singlet = self.peakArr[i].is_singlet
                coster_kronig = self.peakArr[i].is_coster_kronig
                if coster_kronig == True:
                    if self.peakArr[i].peakType.lower() == "double lorentzian":
                        asym_CK = round(self.peakArr[i].asym_CK, 7)
                    else:
                        asym_CK = self.peakArr[i].asymmetry
                    width_CK = self.peakArr[i].lorentz_CK
                else:
                    asym_CK = self.peakArr[i].asymmetry
                    width_CK = self.peakArr[i].lorentz
                branch = self.peakArr[i].branching_ratio
                split = self.peakArr[i].spinOrbitSplit

                FWHM[i] = self.peakArr[i].getFWHM(x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split)

            return FWHM

    def uniqueness(self,x,y, peakType, backgrounds, is_coster_kronig):
            #A function for calculating the uniqueness value of each parameter in a fit. 
            #Outputs plots as txt files for the user to plot elsewhere (too many plots for the program to do?)
            #Still need to get uniqueness values for background parameters...
            self.csv_generate_from = pathlib.Path()
            self.folder_name = pathlib.Path(filedialog.askdirectory(initialdir=pathlib.Path.cwd(), title="Choose a folder"))
            yFit_2 = [0]*len(x)
        
            #Need the total number of parameters for the reduced chi2 values later
            peak_params = 0
            for i in range(self.nPeaks):

                #Getting all parameters
                peakType = self.peakArr[i].peakType
                if(peakType == "Doniach-Sunjic"):
                    asymD = self.peakArr[i].asymmetryDoniach
                else:
                    asymD = 0
                BE = self.peakArr[i].bindingEnergy
                width = self.peakArr[i].lorentz
                sigma = self.peakArr[i].gaussian
                A = self.peakArr[i].amp
                asym = self.peakArr[i].asymmetry
                singlet = self.peakArr[i].is_singlet
                if len(is_coster_kronig) == 0: #Condition for if not coster-kronig
                    for i in range(self.nPeaks):
                        is_coster_kronig.append(False)
                coster_kronig = is_coster_kronig[i]  #List index out of range? During fitting but plots at the end....? Needs to be at least 4 generations done before plotting. peak num = gen num dependent?
                if coster_kronig == True:
                    if self.peakArr[i].peakType.lower() == "double lorentzian":
                        asym_CK = round(self.peakArr[i].asym_CK, 7)
                    else:
                        asym_CK = self.peakArr[i].asymmetry
                    width_CK = self.peakArr[i].lorentz_CK
                
                else:
                    width_CK = self.peakArr[i].lorentz
                    asym_CK = self.peakArr[i].asymmetry
                branch = self.peakArr[i].branching_ratio
                split = self.peakArr[i].spinOrbitSplit




                true_peakType = peakType 
                true_BE = BE
                true_width = width
                true_sigma = sigma
                true_A = A 
                true_asym = asym
                true_singlet = singlet
                true_ck = coster_kronig
                true_Br = branch
                true_SOS = split
                true_asym = asymD
                true_width_CK = width_CK
                true_asym_CK = asym_CK 




                P = 4 # Number of PEAK parameters
                if peakType.lower() == 'double lorentzian' or peakType.lower() == 'doniach-sunjic':
                        P = 5

                if singlet == False:
                    P = 6
                    if peakType.lower() == 'double lorentzian':
                        P = 7
                        if coster_kronig == True:
                            P = 9
                peak_params += P #total number of PEAK parameters

            bkgn_params = 0
            for i in range(self.nBackgrounds):
                B = 1
                if self.bkgn_names == '2-Param Tougaard' or self.bkgn_names == 'Exponential' or self.bkgn_names == 'Polynomial 2':
                    B = 2
                elif self.bkgn_names == '3-Param Tougaard' or self.bkgn_names == 'Polynomial 3':  
                    B = 3
                elif self.bkgn_names == 'Shirley':
                    B = 2 #self.nPeaks 
                bkgn_params += B


            total_params = peak_params + bkgn_params #Number of TOTAL parameters


            #Copies code from get_fit but then does more --> dont want to do in this func because it would take longer for plot function. Have as an Experts tab option?
            for i in range(self.nPeaks):
            
                peakType = self.peakArr[i].peakType
                if(peakType == "Doniach-Sunjic"):
                    asymD = self.peakArr[i].asymmetryDoniach
                else:
                    asymD = 0
                BE = self.peakArr[i].bindingEnergy
                width = self.peakArr[i].lorentz
                sigma = self.peakArr[i].gaussian
                A = self.peakArr[i].amp
                asym = self.peakArr[i].asymmetry
                singlet = self.peakArr[i].is_singlet
                if len(is_coster_kronig) == 0: #Condition for if not coster-kronig
                    for i in range(self.nPeaks):
                        is_coster_kronig.append(False)
                coster_kronig = is_coster_kronig[i]  #List index out of range? During fitting but plots at the end....? Needs to be at least 4 generations done before plotting. peak num = gen num dependent?
                if coster_kronig == True:
                    if self.peakArr[i].peakType.lower() == "double lorentzian":
                        asym_CK = round(self.peakArr[i].asym_CK, 7)
                    else:
                        asym_CK = self.peakArr[i].asymmetry
                    width_CK = self.peakArr[i].lorentz_CK
                
                else:
                    width_CK = self.peakArr[i].lorentz
                    asym_CK = self.peakArr[i].asymmetry
                branch = self.peakArr[i].branching_ratio
                split = self.peakArr[i].spinOrbitSplit
            


                yFit_2 += self.peakArr[i].getY(x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split)

           

            
            
    
            
                print("Peak", i+1)
                unique_params = [BE, width, sigma, A, asym, asymD, width_CK, asym_CK, branch, split] #All the peak parameters we want to vary for uniqueness
                params_range = 1
                param_step = 0.1
           
                j = 0
                while j < (len(unique_params) -1):
              
                    parameter = ''




                    if unique_params[j] == None:
                        j += 1


                    if j == 0: #BE range
                        params_range =  2
                        param_step = 0.01
                        parameter = "BE"
 
                    elif j == 1: #LFWHM range
                        params_range = 0.5
                        param_step = 0.001
                        parameter = "LFWHM"
                  
                    elif j == 2: #GFWHM range
                        params_range =  2
                        param_step = 0.001
                        parameter = "GFWHM"
                  
                    elif j == 3: #Amp range
                        params_range = unique_params[j]*0.2 #10% of amp
                        param_step = unique_params[j]*0.01
                        parameter = "Amplitude"
                  
                    elif j == 4: #asym range
                        params_range = 20
                        param_step = 0.1
                        parameter = "Asymmetry"
                 
                    elif j == 5: #asymD --> never really use...
                        params_range = 3
                        param_step = 0.01
                        parameter = "Asymmetry Doniach-Sunjic"
               
                    elif j ==6: #LFWHM_CK
                        params_range = 1
                        param_step = 0.001
                        parameter = "LFWHM_CK"
                 
                    elif j == 7: #Br range
                        params_range = 0.2
                        param_step = 0.01
                        parameter = "Branching Ratio"
                  
                    elif j == 8: #SOS range
                        params_range = 0.5
                        param_step = 0.01
                        parameter = "Split-Orbit Splitting"
                

                    out_name = "Peak" + str(i+1) + "_" + str(parameter) + "_" + "uniqueness_val" + '.txt'
                    os.chdir(pathlib.Path.cwd().parent)  
                    unique_peak_folder_path = pathlib.Path(self.folder_name)
                    unique_peak_file = unique_peak_folder_path.joinpath(out_name)
                    os.chdir(pathlib.Path.cwd().joinpath('gui'))



               
                    unique_peak_file = open(unique_peak_file, "w")

                    unique_peak_file.write(str(parameter) + "_" + "Peak" + str(i+1))
                    unique_peak_file.write(str("  "))
                    unique_peak_file.write(str("Chi2_r") + "_" + str(parameter) + "_" + "Peak" + str(i+1))
                    unique_peak_file.write(str(" \n"))
                   

                
                    #print("j", j, "val", params[j])
                    count = 0
                    upper_range =  unique_params[j] + params_range
                    lower_range = unique_params[j]-params_range
                    if lower_range < 0:
                        lower_range = 0
                    number_of_points = (upper_range - lower_range)/param_step
                    number_of_points = round(number_of_points)
                    total_chi2_r = [0]*round(number_of_points)
               
                    for unique_params[j] in np.arange(lower_range, upper_range, param_step):

                        if j == 0:
                            BE = unique_params[j]

                        elif j == 1:
                            width =unique_params[j]

                        elif j == 2:
                            sigma = unique_params[j]

                        elif j == 3:
                            A = unique_params[j]

                        elif j == 4:
                            asym = unique_params[j]

                        elif j == 5:
                            asymD = unique_params[j]

                        elif j == 6:
                            width_CK = unique_params[j]

                        elif j == 7:
                            branch = unique_params[j]

                        elif j == 8:
                            split = unique_params[j]

                        chi2 = 0
                        chi2_r = 0
                        yFit_unique = self.peakArr[i].getY(x, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split)
                    
                        for k in range(len(y)):
                            difference = y[k] - yFit_unique[k]
                            chi2 += ((difference**2)/y[k])
                        
                            chi2_r += chi2/(len(y) - total_params)
                        total_chi2_r[j] = (chi2_r)
                     
                        
                        #Figure out way to output data in file. Also run to make sure it work!

                        unique_peak_file.write(str(unique_params[j]))
                        unique_peak_file.write(str("  "))
                        unique_peak_file.write(str(total_chi2_r[j]))
                        unique_peak_file.write(str(" \n"))


                        #print(params[j], total_chi2_r[j])
                    #Save data as (yFit_2, chi2_r)
                    unique_peak_file.close()
                    print(parameter, "uniqueness value saved")

                    #Redefine parameters
                    peakType = true_peakType
                    BE = true_BE
                    width =true_width
                    sigma =true_sigma
                    A = true_A
                    asym = true_asym
                    singlet = true_singlet
                    coster_kronig = true_ck
                    branch = true_Br
                    split = true_SOS
                    asymD = true_asym
                    width_CK = true_width_CK
                    asym_CK = true_asym_CK
        
                    j += 1

          
      







            print("---------------------------")
            print("Output Complete")

            self.tot_peak_fit = yFit_2

            n = 0
            shirley = []*self.nPeaks
            solo_shirley = []
            for l in range(self.nBackgrounds):
                peakType = self.peakArr[n].peakType
                if(peakType == "Doniach-Sunjic"):
                    asymD = self.peakArr[n].asymmetryDoniach
                else:
                    asymD = 0
                BE = self.peakArr[n].bindingEnergy
                width = self.peakArr[n].lorentz
                sigma = self.peakArr[n].gaussian
                A = self.peakArr[n].amp
                asym = self.peakArr[n].asymmetry
                singlet = self.peakArr[n].is_singlet
                coster_kronig = is_coster_kronig[n] #self.peakArr[n].is_coster_kronig
                if coster_kronig == True:
                    if self.peakArr[n].peakType.lower() == "double lorentzian":
                        asym_CK = round(self.peakArr[n].asym_CK, 7)
                    else:
                        asym_CK = self.peakArr[n].asymmetry
                    width_CK = self.peakArr[n].lorentz_CK
                else:
                    asym_CK = self.peakArr[n].asymmetry
                    width_CK = self.peakArr[n].lorentz
                branch = self.peakArr[n].branching_ratio
                split = self.peakArr[n].spinOrbitSplit



                #yFit += self.bkgnArr[i].getY(x,y, BE , width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, self.scale_var) #ISSUE COMES FROM HERE --> NOT SEEING THE self.bindingEnergy or self.lorentz VALUES

                if self.bkgn_names[l] == 'SVSC':
                    self.count_bkgn += 1
                    shirley.append(self.bkgnArr[l].getY(x,y, BE , width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, self.scale_var, peakType, self.tot_peak_fit))
                    n+=1 #Iterate through the various peaks
                    if n == self.nPeaks:  
                        new_shirley = [0]*len(x)    
                        for n in range(self.nPeaks):

                            for i in range(len(x)): #Adding all shirley bkgns together 
                                new_shirley[i] += shirley[n][i] 
                        #for o in range(len(x)):
                        #    new_shirley[o] /= self.nPeaks 
                        yFit_2 += new_shirley
                elif self.bkgn_names[l] == 'Shirley':
                
                    solo_shirley.append(self.bkgnArr[l].getY(x,y, BE , width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, self.scale_var, peakType, self.tot_peak_fit))
                    n+=1 #Iterate through the various peaks
                
                    new_solo_shirley = [0]*len(x)    
      

                    for i in range(len(x)): #Adding all shirley bkgns together 
                        new_solo_shirley[i] += solo_shirley[0][i] 
                
                    yFit_2 += new_solo_shirley

                else:
                    yFit_2 += self.bkgnArr[l].getY(x,y, BE , width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, self.scale_var, peakType, self.tot_peak_fit)
                if self.count_bkgn == self.nPeaks:
                    self.count_bkgn = 0
                    background.get_reset(self,x,y, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, peakType, self.tot_peak_fit)

                elif self.nBackgrounds == self.nPeaks and l + 1 == self.nPeaks: #If only Shirley background
                    background.get_reset(self,x,y, BE, width, sigma, A, asym, asymD, singlet, coster_kronig, width_CK, asym_CK, branch, split, peakType, self.tot_peak_fit)


            #Gives full chi2 values for fit
            #chi2_OG = 0
            #chi2_r_OG = 0
            #for i in range(len(y)):
                #chi2_OG += ((y[i] - yFit_2[i])/y[i])
                #chi2_r_OG += chi2/len(y - total_params) #Reduced chi2 --> chi2/ number of data points - number of parameters in fit

            unique_bkgn_params = [BE, width, sigma, A, asym, asymD, width_CK, asym_CK, branch, split]



            return yFit_2
