# -*- coding: utf-8 -*-
"""
Example setup how to the inverse function to find solutions which fit to the 
existing forward model.

Lotta Ylä-Mella 5.5.2020
"""

from all_functions import inverse
import datetime

import numpy as np


#Parameter values are given as lists to be able to run models with several 
#different limitations

erosion_i = [None, 2]       # Erosion in the inverse solution [m]
tie = [None, 10]            # Tied last deglaciation [ka]
mini_ex = [1,2]             # Minimum number of exposures
maxi_ex = [4,2]             # Maximum number of exposures
rand_maxi = [50,50]         # Maximum time of the models
maxi_ero =[4,10]            # Maximum erosion
depth = [0,2]               # Sample depth
names=['Test1', 'Test2']    # Filenames


for k in range(len(names)):
    
    # Used isotope
    isotope = 1
    
    # Starting times of glaciations (ice) and deglaciations in the 
    # original forward model; reference model
    time_ice_fwd = np.array([25,0])     
    time_degla_fwd = np.array([35,10])  
    
    # Erosion in the forward model
    erosion_fwd = np.array([2])
    
    # Tied last deglaciation
    tied=tie[k]
    
    # Minimum and maximum number of exposures
    min_ex=mini_ex[k]
    max_ex=maxi_ex[k]
    
    # Maximum time of models
    rand_max = rand_maxi[k]
    
    # Maximum erosion
    max_erosion = maxi_ero[k]
    
    # Sample depths
    z_chosen = np.array([depth[k]])
    
    # Sample error
    z_error = 0.025
    
    # Erosion in the inverse models
    erosion_inv = np.array([erosion_i[k]])
    
    # Possibility to preset times for inversion
    time_ice_inv = None
    time_degla_inv = None
    
    # Maximum numbers of exposures set to 4 to make files to have similar 
    # structure even if they would have different number of exposures
    max_complexity = 4
    
    #Cols for times*2, erosion (one less), isotope, complexity and misfit
    columns = max_complexity*2 + max_complexity-1 + 3 
    
    # Number of tested solutions
    n = 10
    
    # Array to save solutions
    models = np.zeros((n,columns))
    
    
    
    
    for i in range(n):
        # Randomly choose how complicated model is tested
        complexity = np.random.randint(min_ex,max_ex+1)

        # Find inverse solution
        free_para, all_info, misfit, misfit_arr, chunk_erosion = \
        inverse(z_chosen, z_error, isotope, time_ice_fwd, time_degla_fwd, 
                erosion_fwd, time_ice_inv, time_degla_inv, erosion_inv, 
                complexity, tied, rand_max, max_erosion)

        # Fill models array to be indepent of complexity of the solution
        # Add times to the array
        for j in range(max_complexity):
            if complexity == 4:
                models[i,j] = all_info[j]
                models[i,j+max_complexity] = all_info[j+complexity]
            elif complexity == 3:
                if j < 1:
                    models[i,j] = 0
                    models[i,j+max_complexity] = 0
                else:
                    models[i,j] = all_info[j-1]
                    models[i,j+max_complexity] = all_info[j+complexity-1]
            elif complexity == 2:
                if j < 2:
                    models[i,j] = 0
                    models[i,j+max_complexity] = 0
                else:
                    models[i,j] = all_info[j-2]
                    models[i,j+max_complexity] = all_info[j+complexity-2]
            elif complexity == 1:
                if j < 3:
                    models[i,j] = 0
                    models[i,j+max_complexity] = 0
                else:
                    models[i,j] = all_info[j-3]
                    models[i,j+max_complexity] = all_info[j+complexity-3]
        
        # Add erosion to the array
        for j in range(len(chunk_erosion)):
            models[i,j+2*max_complexity] = chunk_erosion[j]#erosion_inv[j]
          
        # Add other parameters
        models[i,-3] = free_para[0] #Isotope
        models[i,-2] = complexity
        models[i,-1] = misfit
        
        # If needed, record the time of how long the code runs
        #if i%200 == 0:
        #    print(i,datetime.datetime.now())
    

    # Save the results to a file
    # File name
    name = names[k]
    
    # Set all the variables in the beginning of the file
    setup = {'Isotope':isotope, 'Start of ice':time_ice_fwd,
             'Start of exp':time_degla_fwd, 'Erosion':erosion_fwd,
             'Sample depths':z_chosen, 'Error in depth':z_error,
             'Maximum exposures':max_ex, 'Tied':tied, 
             'Inverse erosion':erosion_inv, 'Total maximum time':rand_max, 
             'Maximum erosion':max_erosion}
   
    # Save the results with a header
    np.savetxt(name+'.txt',models,header=str(setup))
    
