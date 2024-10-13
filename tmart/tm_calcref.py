# This file is part of T-Mart.
#
# Copyright 2023 Yulun Wu.
#
# T-Mart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import numpy as np 
import pandas as pd
import sys
from copy import copy

# Analyze the output of TMart and differentiate direct, env and atm intrinsic reflectances  
def calc_ref(df, n_photon = None, detail = False):
    '''Analyze the results of T-Mart and calculate reflectances. 
    
    Arguments:

    * ``df`` -- Results from T-Mart runs
    * ``n_photon`` -- Specify the number of photons in the run when firing the photon upwards. If not specified, the number of unique pt_id will be used. This can lead to errors when photons were fired upwards because some photons will not have pt_id.
    * ``detail`` -- Boolean. Differentiate Cox-Munk, whitecap, water-leaving and land contributions

    Output:

    * A list of atmospheric intrinsic reflectance, direct reflectance, environmental reflectance and total reflectance. 

    Example usage::

      R = tmart.calc_ref(results)
      for k, v in R.items():
          print(k, '\t ' , v)

    '''
    
    print('=====================================')
    print('Calculating radiometric quantities...')
    
    # Columes: 0 pt_id, 1 movement, 2 L_cox-munk, 3 L_whitecap, 4 L_water, 5 L_land, 
    # 6 L_rayleigh, 7 L_mie, 8 9 10 surface xyz, 11 shadowed, 12 if_env
    
    if n_photon == None:
        n_photon = np.unique(df[:,0]).shape[0]
    
    R_atm = np.sum(df[:,6:8]) / n_photon
    R_dir = np.sum(df[df[:,12] == 0 ,2:6]) / n_photon # if_env == 1 and all surface reflectance 
    R_env = np.sum(df[df[:,12] == 1 ,2:6]) / n_photon # if_env == 0 and all surface reflectance 
    R_total = R_atm + R_dir + R_env
    
    if detail:
        R_dir_coxmunk   = np.sum(df[df[:,12] == 0 ,2]) / n_photon
        R_dir_whitecap  = np.sum(df[df[:,12] == 0 ,3]) / n_photon
        R_dir_water     = np.sum(df[df[:,12] == 0 ,4]) / n_photon
        R_dir_land      = np.sum(df[df[:,12] == 0 ,5]) / n_photon
        
        R_env_coxmunk   = np.sum(df[df[:,12] == 1 ,2]) / n_photon
        R_env_whitecap  = np.sum(df[df[:,12] == 1 ,3]) / n_photon
        R_env_water     = np.sum(df[df[:,12] == 1 ,4]) / n_photon
        R_env_land      = np.sum(df[df[:,12] == 1 ,5]) / n_photon        
        
        R_output = {'R_atm':R_atm,
                    'R_dir':R_dir,
                    '_R_dir_coxmunk':R_dir_coxmunk,
                    '_R_dir_whitecap':R_dir_whitecap,
                    '_R_dir_water':R_dir_water,
                    '_R_dir_land':R_dir_land,
                    'R_env':R_env,
                    '_R_env_coxmunk':R_env_coxmunk,
                    '_R_env_whitecap':R_env_whitecap,
                    '_R_env_water':R_env_water,
                    '_R_env_land':R_env_land,
                    'R_total':R_total}
    
    else:      
        R_output = {'R_atm':R_atm,
                    'R_dir':R_dir,
                    'R_env':R_env,
                    'R_total':R_total}
    
    
    return R_output
