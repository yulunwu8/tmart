# This file is part of T-Mart.
#
# Copyright 2023 Yulun Wu.
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

# Specify wavelength in nm
wl = 800

# Distance to shore in m 
dist = 1000

# Built-in spectral libraries 
water = tmart.SpectralSurface('water_chl1')
vegetation = tmart.SpectralSurface('vegetation')


### DEM and reflectance ###

# Three same-size numpy arrays are needed
image_DEM = np.full((2, 2), 0) # in meters
image_reflectance = np.array([[vegetation.wl(wl),vegetation.wl(wl)],
                              [water.wl(wl),water.wl(wl)]]) 
image_isWater = np.array([[0,0],[1,1]])

# pixel width and length, in meters 
cell_size = 20_000

# Synthesize a surface object
my_surface = tmart.Surface(image_DEM,image_reflectance,image_isWater,cell_size)  
# Set background information, 1 or 2 background surfaces can be set;
# If 2: the first background is the one closer to [0,0]
my_surface.set_background(bg_ref        = [vegetation.wl(wl),water.wl(wl)], 
                          bg_isWater    = [0,1], # water
                          bg_elevation  = 0, # elevation of both background
                          bg_coords     = [[0,20_000],[40_000,20_000]]) # a line dividing two background                                    

### Atmosphere ###

# Atmophere profile comes from 6S
atm_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeSummer) 
aerosol_type = 'Maritime' 
aot550 = 0.1

# Synthesize an atmosphere object    
my_atm = tmart.Atmosphere(atm_profile, aot550, aerosol_type)   

### Running T-Mart ###

# Make a T-Mart object 
my_tmart = tmart.Tmart(Surface = my_surface, Atmosphere= my_atm, shadow=False)
my_tmart.set_wind(wind_speed=10, wind_azi_avg = True)
my_tmart.set_water(water_salinity=35, water_temperature=20)

# Specify the sensor's position (x, y, z), viewing direction relative 
# to the sensor (zenith, azimuth), sun's direction relative to the target 
# (zenith, azimuth)
my_tmart.set_geometry(target_coords=[30000,20000+dist], 
                      target_pt_direction=[150,0],
                      sun_dir=[30,90])  

# Run single photon and plot 
# results = my_tmart.run_plot(wl=wl, plot_on=True, plot_range=[0,100_000,0,100_000,0,100_000])

if __name__ == "__main__":
    
    n_photon = 100_000
    results = my_tmart.run(wl=wl, n_photon=n_photon)
    
    # Calculate reflectances using recorded photon information 
    R = tmart.calc_ref(results, detail=True)
    
    for k, v in R.items():
        print(k, '     ' , v)
    
    
    
    
    