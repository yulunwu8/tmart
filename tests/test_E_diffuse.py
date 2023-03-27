# This file is part of TMart.
#
# Copyright 2022 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# %matplotlib qt

import sys
sys.path.append('/Users/yw/Desktop/tmart')

import tmart
import numpy as np
import pandas as pd
from Py6S.Params.atmosprofile import AtmosProfile




def calculate_Ediffuse(albedo, ot_abs, ot_rayleigh, sza, n_photon):
    
    # Specify wavelength in nm, not used as atm is specified
    wl = 400
    
    ### DEM and reflectance ###
    
    # Three same-size numpy arrays are needed
    image_DEM = np.full((2, 2), 0) # in meters
    image_reflectance = np.full((2, 2), albedo) # unitless     
    image_isWater = np.full((2, 2), 0) # 1 is water, 0 is land
    
    # pixel width and length, in meters 
    cell_size = 20_000 
    
    # Synthesize a surface object
    my_surface = tmart.Surface(DEM = image_DEM,
                               reflectance = image_reflectance,
                               isWater = image_isWater,
                               cell_size = cell_size)  
    
    ### Atmosphere ###
    
    # Atmophere profile comes from 6S
    atm_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeSummer) 
    
    # Synthesize an atmosphere object 
    my_atm = tmart.Atmosphere(atm_profile,n_layers=1,
                              specify_abs=ot_abs,
                              specify_ot_rayleigh=ot_rayleigh)
    
    ### Running T-Mart ###
    
    # Make a T-Mart object 
    my_tmart = tmart.Tmart(Surface = my_surface, Atmosphere= my_atm)
    
    # Specify the sensor's position (x, y, z), viewing direction relative 
    # to the sensor (zenith, azimuth), sun's direction relative to the target 
    # (zenith, azimuth), see Parameters for more geometry options
    my_tmart.set_geometry(sensor_coords=[51,50,0.00001], 
                          target_pt_direction='lambertian_up',
                          sun_dir=[sza,0])
    
    results = my_tmart.run(wl=wl,n_photon=n_photon,njobs=2000)
    
    # Calculate reflectances using recorded photon information 
    R = tmart.calc_ref(results, n_photon=n_photon)
    
    print('R: '+str(R))


    return R
    




albedo = 0.1

ot_abs = 0.0

n_photon = 10_000


ot_rayleigh = 0.36
sza = 30
R = calculate_Ediffuse(albedo, ot_abs, ot_rayleigh, sza, n_photon)




'''
results = []

for sza in range(0,100,20):
    
    for ot_rayleigh in np.arange(0.05, 0.55, 0.05):
        print ('\n------------------------------------ SZA: ' + str(sza))
        print ('------------------------------------ ot_rayleigh: ' + str(ot_rayleigh))
        
        R = calculate_Ediffuse(albedo, ot_abs, ot_rayleigh, sza, n_photon)
        results.append({'sza': sza, 'ot_rayleigh': ot_rayleigh, 'E_diffuse': R})

results = pd.DataFrame(results)



results.to_csv('E_diffuse.csv', index=False)

'''




























