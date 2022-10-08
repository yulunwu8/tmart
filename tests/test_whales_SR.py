# This file is part of TMart.
#
# Copyright 2022 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.








# To calculate the glint reflectance in whale images 

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


# 460 530 590

# Band 2 of Sentinel-2A and its central wavelength
band = Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_02)
band = None
wl = 590 

### DEM and reflectance ###

# Three same-size numpy arrays are needed
image_DEM = np.full((2, 2), 0) # in meters
image_reflectance = np.full((2, 2), 0) # unitless     
image_isWater = np.full((2, 2), 1) # 1 is water, 0 is land

# pixel width and length, in meters 
cell_size = 20_000 

# Synthesize a surface object
my_surface = tmart.Surface(DEM = image_DEM,
                           reflectance = image_reflectance,
                           isWater = image_isWater,
                           cell_size = cell_size)  
# Set background information, 1 or 2 background surfaces can be set;
# If 2: the first background is the one closer to [0,0]
my_surface.set_background(bg_ref        = 0, # background reflectance
                          bg_isWater    = [1,1], # if is water
                          bg_elevation  = 0, # elevation of both background
                          bg_coords     = [[0,0],[10,10]]) # a line dividing the two background                                    

### Atmosphere ###

# Atmophere profile comes from 6S
atm_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeSummer) 
aerosol_type = 'Maritime' 
aot550 = 0.1
n_layers = 20
aerosol_scale_height = 2 # Unless you have a reason, don't change this

# Synthesize an atmosphere object    
my_atm = tmart.Atmosphere(atm_profile, aot550, aerosol_type, n_layers, aerosol_scale_height)

### Running T-Mart ###

# Make a T-Mart object 
my_tmart = tmart.Tmart(Surface = my_surface, Atmosphere= my_atm, shadow=True)
my_tmart.set_wind(wind_speed=10, wind_azi_avg = True)
my_tmart.set_water(water_salinity=35, water_temperature=20)

# Specify the sensor's position (x, y, z), viewing direction relative 
# to the sensor (zenith, azimuth), sun's direction relative to the target 
# (zenith, azimuth), see Parameters for more geometry options
my_tmart.set_geometry(sensor_coords=[51,50,1], 
                      target_pt_direction=[150,0],
                      sun_dir=[30,180])

# Number of photons
n_photon = 10_000

results = my_tmart.run(wl=wl, band=band, n_photon=n_photon)
# Calculate reflectances using recorded photon information 
R = tmart.calc_ref(results, detail=True)
for k, v in R.items():
    print(k, '     ' , v)

















