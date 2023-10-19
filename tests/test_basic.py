# This file is part of TMart.
#
# Copyright 2022 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# %matplotlib qt


# Go up by 2 directory and import 

import sys
import os.path as path
two_up =  path.abspath(path.join(__file__ ,"../.."))
sys.path.append(two_up)


import tmart
import numpy as np
from Py6S.Params.atmosprofile import AtmosProfile
import Py6S

# Specify wavelength in nm
band = Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_08)
wl = 833

### DEM and reflectance ###
image_DEM = np.array([[0,0],[0,0]]) # in meters
image_reflectance = np.array([[0.00,0.00],[0.00,0.00]]) # unitless     
# image_isWater = np.zeros(image_DEM.shape) # 1 is water, 0 is land
image_isWater = np.array([[1,1],[1,1]]) 


# Synthesize a surface object
my_surface = tmart.Surface(DEM = image_DEM,
                           reflectance = image_reflectance,
                           isWater = image_isWater,
                           cell_size = 10_000)  
           
my_surface.set_background(bg_ref        = 0.0, # background reflectance
                          bg_isWater    = 1, # if is water
                          bg_elevation  = 0, # elevation of both background
                          bg_coords     = [[0,0],[10,10]]) # a line dividing the two background                                    
                    
### Atmosphere ###
atm_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeSummer) 


my_atm = tmart.Atmosphere(atm_profile, aot550 = 0.11275386685706923, aerosol_type = 0.287 )

### Running T-Mart ###
my_tmart = tmart.Tmart(Surface = my_surface, Atmosphere= my_atm, shadow=False)
my_tmart.set_wind(wind_speed=1, wind_azi_avg = True)

sensor_coords=[51,50,130_000]


my_tmart.set_geometry(sensor_coords=sensor_coords, 
                      target_pt_direction=[170.52804413432926, 191.91873559828522],
                      sun_dir=[30.9608405674786, 323.9885587375248])

# my_tmart.set_geometry(sensor_coords=sensor_coords, 
#                       target_pt_direction=[180-30, 0],
#                       sun_dir=[10, 0])


results = my_tmart.run(wl=wl, band=band, n_photon=10_000)

# Calculate reflectances using recorded photon information 
R = tmart.calc_ref(results,detail=True)
for k, v in R.items():
    print(k, '     ' , v)



my_atm_profile = my_tmart.atm_profile_wl
np.sum(my_atm_profile.ot_mie)











