from PhysicsModules.XPS.xps_neo.helper import str_to_bool
# ------------------------------------
# Andy Lau
# 4/22/2019
# goal : Parsed the ini into each specific documents.
#-------------------------------------

#UPDATE: ANY CHANGES IN THIS FILE WILL ONLY UPDATE WITH pip install . to the cd XPS_Neo location

def split_string(var_dict,label):
    """
    Read the path list
    """
    # print(num_compounds)
    arr_str = var_dict[label]
    starter = []
    end = []
    k = 0
    split_str = []
    for i in arr_str:
        if i == '[':
            starter.append(k)
        elif i == ']':
            end.append(k)
        k = k + 1


    assert(len(starter) == len(end)),'Bracket setup not right.'
    # if num_compounds > 1:
    # 	assert(num_compounds == len(starter)),'Number of compounds not matched.'
    # 	assert(num_compounds == len(end)),'Number of compounds not matched.'

    # check if both are zeros, therefore the array is one 1 dimensions
    # arr_str = optional_var(var_dict,label,[],list)
    if len(starter) == 0 and len(end) == 0:
        split_str = list(arr_str.split(","))
    f_arr = []
    for i in range(len(split_str)):
        f_arr.append(float(split_str[i]))

    return f_arr

def optional_var(dict,name_var,alt_var=None,type_var=int):
    """
    Detections of optional variables exists within input files, and
        put in corresponding default inputs parameters.
    """
    # boolean needs special attentions
    if type_var == bool:
        if name_var in dict:
            return_var = str_to_bool(dict[name_var])
        else:
            return_var = alt_var
    elif type_var is None:
        if name_var in dict:
            return_var = dict[name_var]
        else:
            return_var = None
    else:
        if name_var in dict:
            return_var = type_var(dict[name_var])
        else:
            return_var = type_var(alt_var)
    return return_var

def optional_range(var_dict,label):
    if label not in var_dict:
        # Label is not therefore
        return_var = []
    else:
        return_var = split_string(var_dict,label)

    return return_var

def optional_range_string(var_dict, label):
    if label not in var_dict:
        return_var = []
    else:
       return_var = var_dict[label].split(',')


    return return_var



def load_config(file_dict):
    """Parse the INI section dict into config values (Phase 2).

    The body below is the former module-level parsing code, unchanged,
    so value computation and print order are identical. The returned
    dict's per-peak lists are shared mutable objects by design -
    XPS_GA.addPeak/removePeak mutate them in place (design.md Phase 2).
    """
    # -----
    Inputs_dict = file_dict['Inputs']
    Populations_dict = file_dict['Populations']
    Mutations_dict = file_dict['Mutations']
    Paths_dict = file_dict['Paths']
    Outputs_dict = file_dict['Outputs']
    # -----

    # Input
    data_file = Inputs_dict['data_file']
    output_file = Inputs_dict['output_file']
    skipLn = int(Inputs_dict['skipln'])
    x_offset = float(Inputs_dict['x_offset'])
    y_offset = float(Inputs_dict['y_offset'])

    #data_cutoff = [float(x) for x in Inputs_dict['data_cutoff'].split(',')] -Unneccesary, remove later when we have stable build
    #pathrange_file = optional_var(Inputs_dict,'pathrange_file',None,None)		-evan


    # population
    size_population = int(Populations_dict['population'])
    number_of_generation = int(Populations_dict['num_gen'])
    best_sample = int(Populations_dict['best_sample'])
    lucky_few = int(Populations_dict['lucky_few'])

    # Mutations
    chance_of_mutation = int(Mutations_dict['chance_of_mutation'])
    original_chance_of_mutation = int(Mutations_dict['original_chance_of_mutation'])
    mutated_options = int(Mutations_dict['mutated_options'])

    # Paths  -Should change the name of this to be peak or params or something later for readability - evan
    npaths = int(Paths_dict['npeaks'])


    #Rearranged to put all input guesses in order they appear in Fitting Parameters tab
    BE_guess = optional_range(Paths_dict,'be')
    sigma_guess = optional_range(Paths_dict,'sigma_guess')
    gamma_guess = optional_range(Paths_dict,'gamma_guess')
    amp_guess = optional_range(Paths_dict,'amp_guess')
    asym_guess = optional_range(Paths_dict,'asym_guess')
    sos_guess = optional_range(Paths_dict,'sos_guess')
    br_guess = optional_range(Paths_dict,'br_guess')




    is_singlet = optional_range_string(Paths_dict,'is_singlet')
    for i in range(len(is_singlet)):
        if is_singlet[i].strip() == 'True':
            is_singlet[i] = True
        else:
            is_singlet[i] = False
    #branching_ratio = optional_range(Paths_dict,'branching_ratio')

    #spinOrbitSplit_guess = optional_range(Paths_dict,'spinorbitsplit') #READS IT ALL IN LOWERCASE --> OTHERWISE YOU GET AN ERROR

    is_coster_kronig = optional_range_string(Paths_dict,'is_coster_kronig')
    for i in range(len(is_coster_kronig)):
        if is_coster_kronig[i].strip() == 'True':
            is_coster_kronig[i] = True
        else:
            is_coster_kronig[i] = False
    alt_lorentz = optional_range_string(Paths_dict,'alt_LFWHM')
    peak_add_remove = optional_range_string(Paths_dict,'peak_adding')
    scale_data = optional_range_string(Paths_dict,'scale_bool')

    for i in range(len(alt_lorentz)):
        if alt_lorentz[i] == 'True':
            alt_lorentz = True
        else:
            alt_lorentz[i] = False

    for i in range(len(peak_add_remove)):
        if peak_add_remove[i] == 'True':
            peak_add_remove = True
        else:
            peak_add_remove[i] = False

    for i in range(len(scale_data)):
        if scale_data[i] == 'True':
            scale_data = True
        else:
            scale_data[i] = False

    gen_alt = int(Paths_dict['gen_alt_val'])




    element = str(Paths_dict['element_select'])
    photoelectronLine = str(Paths_dict['photoline_select'])

    #spinOrbitSplit_guess = optional_range(Paths_dict,'spinorbitsplit')

    #Coster-Kronig Values for higher BE peak:
    #Older / hand-made INIs may omit these keys or give fewer per-peak guesses
    #than there are peaks (Individual indexes them per peak, so a short list
    #raises IndexError). Fall back to the GUI defaults for the range and copy
    #the sibling per-peak guess list when the CK list is missing or too short.
    #The copy is a distinct object because addPeak/removePeak mutate these
    #lists in place, matching xps.py's asym_CK_guess.append(asym_guess[-1]).
    gamma_CK_range = optional_range(Paths_dict,'gamma_ck_range')
    if not gamma_CK_range:
        gamma_CK_range = [-0.01, 2.0, 0.0001]
    gamma_CK_guess = optional_range(Paths_dict,'gamma_ck_guess') #Same as other gamma to start
    if len(gamma_CK_guess) < npaths:
        gamma_CK_guess = list(gamma_guess)

    asym_CK_range = optional_range(Paths_dict,'asym_ck_range')
    if not asym_CK_range:
        asym_CK_range = [1.0, 20.0, 0.001]
    asym_CK_guess = optional_range(Paths_dict,'asym_ck_guess') #Same as other asym to start
    if len(asym_CK_guess) < npaths:
        asym_CK_guess = list(asym_guess)

    #branching_ratio_range = optional_range(Paths_dict,'branching_ratio_range')
    #spinOrbitSplit_range = optional_range(Paths_dict,'spinorbitsplit_range')
    background_type = optional_range_string(Paths_dict,'background_type')
    peak_type = optional_range_string(Paths_dict,'peak_type')
    #BE_range = optional_range(Paths_dict,'be_range')
    BE_range_min = optional_range(Paths_dict,'be_range_min')
    BE_range_max = optional_range(Paths_dict,'be_range_max')
    BE_range_delta = optional_range(Paths_dict,'be_range_delta')
    BE_limited = optional_range_string(Paths_dict,'be_limited')
    for i in range(len(BE_limited)):
        if BE_limited[i].strip() == 'True':
            BE_limited[i] = True
        else:
            BE_limited[i] = False
    BE_correlated = optional_range_string(Paths_dict,'be_correlated')
    for i in range(len(BE_correlated)):
        if BE_correlated[i].strip() == 'Peak #':
            BE_correlated[i] = i + 1
        else:
            BE_correlated[i] = BE_correlated[i]
    BE_correlated_mult = optional_range(Paths_dict,'be_correlated_mult')

    sigma_range_min = optional_range(Paths_dict,'sigma_range_min')
    sigma_range_max = optional_range(Paths_dict,'sigma_range_max')
    sigma_range_delta = optional_range(Paths_dict,'sigma_range_delta')
    sigma_limited = optional_range_string(Paths_dict,'sigma_limited')
    for i in range(len(sigma_limited)):
        if sigma_limited[i].strip() == 'True':
            sigma_limited[i] = True
        else:
            sigma_limited[i] = False
    sigma_correlated = optional_range_string(Paths_dict,'sigma_correlated')
    for i in range(len(sigma_correlated)):
        if sigma_correlated[i].strip() == 'Peak #':
            sigma_correlated[i] = i + 1
        else:
            sigma_correlated[i] = sigma_correlated[i]
    sigma_correlated_mult = optional_range(Paths_dict,'sigma_correlated_mult')

    gamma_range_min = optional_range(Paths_dict,'gamma_range_min')
    gamma_range_max = optional_range(Paths_dict,'gamma_range_max')
    gamma_range_delta = optional_range(Paths_dict,'gamma_range_delta')
    gamma_limited = optional_range_string(Paths_dict,'gamma_limited')
    for i in range(len(gamma_limited)):
        if gamma_limited[i].strip() == 'True':
            gamma_limited[i] = True
        else:
            gamma_limited[i] = False
    gamma_correlated = optional_range_string(Paths_dict,'gamma_correlated')
    for i in range(len(gamma_correlated)):
        if gamma_correlated[i].strip() == 'Peak #':
            gamma_correlated[i] = i + 1
        else:
            gamma_correlated[i] = gamma_correlated[i]
    gamma_correlated_mult = optional_range(Paths_dict,'gamma_correlated_mult')

    amp_range_min = optional_range(Paths_dict,'amp_range_min')
    amp_range_max = optional_range(Paths_dict,'amp_range_max')
    amp_range_delta = optional_range(Paths_dict,'amp_range_delta')
    amp_limited = optional_range_string(Paths_dict,'amp_limited')
    for i in range(len(amp_limited)):
        if amp_limited[i].strip() == 'True':
            amp_limited[i] = True
        else:
            amp_limited[i] = False
    amp_correlated = optional_range_string(Paths_dict,'amp_correlated')
    for i in range(len(amp_correlated)):
        if amp_correlated[i].strip() == 'Peak #':
            amp_correlated[i] = i + 1
        else:
            amp_correlated[i] = amp_correlated[i]
    amp_correlated_mult = optional_range(Paths_dict,'amp_correlated_mult')



    asym_range_min = optional_range(Paths_dict,'asym_range_min')
    asym_range_max = optional_range(Paths_dict,'asym_range_max')
    asym_range_delta = optional_range(Paths_dict,'asym_range_delta')
    asym_limited = optional_range_string(Paths_dict,'asym_limited')
    for i in range(len(asym_limited)):
        if asym_limited[i].strip() == 'True':
            asym_limited[i] = True
        else:
            asym_limited[i] = False
    asym_correlated = optional_range_string(Paths_dict,'asym_correlated')
    for i in range(len(asym_correlated)):
        if asym_correlated[i].strip() == 'Peak #':
            asym_correlated[i] = i + 1
        else:
            asym_correlated[i] = asym_correlated[i]
    asym_correlated_mult = optional_range(Paths_dict,'asym_correlated_mult')


    sos_range_min = optional_range(Paths_dict,'sos_range_min')
    sos_range_max = optional_range(Paths_dict,'sos_range_max')
    sos_range_delta = optional_range(Paths_dict,'sos_range_delta')
    sos_limited = optional_range_string(Paths_dict,'sos_limited')
    for i in range(len(sos_limited)):
        if sos_limited[i].strip() == 'True':
            sos_limited[i] = True
        else:
            sos_limited[i] = False
    sos_correlated = optional_range_string(Paths_dict,'sos_correlated')
    for i in range(len(sos_correlated)):
        if sos_correlated[i].strip() == 'Peak #':
            sos_correlated[i] = i + 1
        else:
            sos_correlated[i] = sos_correlated[i]
    sos_correlated_mult = optional_range(Paths_dict,'sos_correlated_mult')




    br_range_min = optional_range(Paths_dict,'br_range_min')
    br_range_max = optional_range(Paths_dict,'br_range_max')
    br_range_delta = optional_range(Paths_dict,'br_range_delta')
    br_limited = optional_range_string(Paths_dict,'br_limited')
    for i in range(len(br_limited)):
        if br_limited[i].strip() == 'True':
            br_limited[i] = True
        else:
            br_limited[i] = False
    br_correlated = optional_range_string(Paths_dict,'br_correlated')
    for i in range(len(br_correlated)):
        if br_correlated[i].strip() == 'Peak #':
            br_correlated[i] = i + 1
        else:
            br_correlated[i] = br_correlated[i]
    br_correlated_mult = optional_range(Paths_dict,'br_correlated_mult')





    #asymmetry_range = optional_range(Paths_dict, 'asymmetry_range')
    asymmetryDoniach_range = optional_range(Paths_dict, 'asymmetrydoniach_range')
    k_range = optional_range(Paths_dict,'k_range')
    background_range = optional_range(Paths_dict,'background_range',)
    #CTou3_range = optional_range(Paths_dict,'ctou3_range')
    #DTou3_range = optional_range(Paths_dict,'ctou3_range')
    #background_shir_range = optional_range(Paths_dict,'background_shir_range',)
    slope_range = optional_range(Paths_dict,'slope_range',)
    #exp_amp_range = optional_range(Paths_dict,'exp_amp_range',)
    #exp_decay_range = optional_range(Paths_dict,'exp_decay_range',)

    for i in range(len(peak_type)):
        peak_type[i].replace(" ","")

    print("Element selected is" +  ' ' + str(element))
    print("Photoelectron Line selected is" + ' ' + str(photoelectronLine))
    print("Background type is " + ', '.join(background_type))
    print("Peak Type is " + ', '.join(peak_type))
    print("Is a Singlet: " + str(is_singlet))
    print("Is a Coster-Kronig: " + str(is_coster_kronig))

    #print("Branching Ratio Range is " + str(branching_ratio_range))
    #print("Branching Ratio is " + str(branching_ratio))
    print("Branching Ratio Limited " + str(br_limited))
    print("Branching Ratio Correlated Peak " + str(br_correlated))
    print("Branching Ratio Correlated Multipliear " + str(br_correlated_mult))
    print("Branching Ratio Range Minimum is " +  str(br_range_min))
    print("Branching Ratio Range Maximum is " +  str(br_range_max))
    print("Branching Ratio Range Delta is " +  str(br_range_delta))
    print("Trying Branching Ratio " + str(br_guess))
    #print("Spin-Orbit Splitting Range is " + str(spinOrbitSplit_range))
    #print("Trying Spin-Orbit Splitting " + str(spinOrbitSplit_guess))
    #print("Binding Energy Range is " + str(BE_range))
    print("Binding Energy Limited " + str(BE_limited))
    print("Binding Energy Correlated Peak " + str(BE_correlated))
    print("Binding Energy Correlated Multipliear " + str(BE_correlated_mult))
    print("Binding Energy Range Minimum is " +  str(BE_range_min))
    print("Binding Energy Range Maximum is " +  str(BE_range_max))
    print("Binding Energy Range Delta is " +  str(BE_range_delta))
    print("Trying BE " + str(BE_guess))
    print("Sigma Limited " + str(sigma_limited))
    print("Sigma Correlated Peak" + str(sigma_correlated))
    print("Sigma Correlated Multipliear " + str(sigma_correlated_mult))
    print("Sigma Range Minimum is " + str(sigma_range_min))
    print("Sigma Range Maximum is " +  str(sigma_range_max))
    print("Sigma Range Delta is " +  str(sigma_range_delta))
    print("Trying Sigma " + str(sigma_guess))
    print("Gamma Limited " + str(gamma_limited))
    print("Gamma Correlated Peak " + str(gamma_correlated))
    print("Gamma Correlated Multipliear " + str(gamma_correlated_mult))
    print("Gamma Range Minimum is " + str(gamma_range_min))
    print("Gamma Range Maximum is " +  str(gamma_range_max))
    print("Gamma Range Delta is " +  str(gamma_range_delta))
    print("Trying Gamma " + str(gamma_guess))


    #IDK if we should state these. They are the same as gamma
    #print("Gamma Coster-Kronig Range is " + str(gamma_CK_range))
    #print("Trying Gamma " + str(gamma_guess))
    print("Amplitude Limited " + str(amp_limited))
    print("Amplitude Correlated Peak " + str(amp_correlated))
    print("Amplitude Correlated Multipliear " + str(amp_correlated_mult))
    print("Amplitude Range Minimum is " +  str(amp_range_min))
    print("Amplitude Range Maximum is " +  str(amp_range_max))
    print("Amplitude Range Delta is " +  str(amp_range_delta))
    print("Trying Amplitude " + str(amp_guess))

    print("Asymmetry Limited " + str(asym_limited))
    print("Asymmetry Correlated Peak " + str(asym_correlated))
    print("Asymmetry Correlated Multipliear " + str(asym_correlated_mult))
    print("Asymmetry Range Minimum is " +  str(asym_range_min))
    print("Asymmetry Range Maximum is " +  str(asym_range_max))
    print("Asymmetry Range Delta is " +  str(asym_range_delta))
    print("Trying Asymmetry " + str(asym_guess))

    print("Spin-Orbit Splitting Limited " + str(sos_limited))
    print("Spin-Orbit Splitting Correlated Peak " + str(sos_correlated))
    print("Spin-Orbit Splitting Correlated Multipliear " + str(sos_correlated_mult))
    print("Spin-Orbit Splitting Range Minimum is " +  str(sos_range_min))
    print("Spin-Orbit Splitting Range Maximum is " +  str(sos_range_max))
    print("Spin-Orbit Splitting Range Delta is " +  str(sos_range_delta))
    print("Trying Spin-Orbit Splitting " + str(sos_guess))

    #print("Asymmetry Range is " + str(asymmetry_range))
    print("Asymmetry Doniach-Sunjic Range is " + str(asymmetryDoniach_range))
    print("K range is " + str(k_range))
    print("Background Range is " + str(background_range))
    #print("C Range is " + str(CTou3_range))
    #print("D Range is " + str(DTou3_range))
    #print("Shirley Background Range is " + str(background_shir_range))
    print("Slope Range is " + str(slope_range))
    #print("Exponential Amplitude Range is " + str(exp_amp_range))
    #print("Exponential Decay Rate Range is " + str(exp_decay_range))
    #print("Generation for altering chi^2: " + str(alt_gen))
    print("ALLOWING LFWHM Alteration: " + str(alt_lorentz))
    print("ALLOWING PEAK ADDITION/SUBTRACTION: " + str(peak_add_remove))
    print("SCALE DATA: " + str(scale_data))
    # Output
    printgraph = str_to_bool(Outputs_dict['print_graph'])
    num_output_paths = str_to_bool(Outputs_dict['num_output_paths']) #PYVAR --> Fix
    steady_state = optional_var(Outputs_dict,'steady_state_exit',False,bool) #PYVAR --> Fix

    config = dict(locals())
    del config['file_dict']
    return config
