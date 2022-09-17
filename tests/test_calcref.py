# This file is part of TMart.
#
# Copyright 2022 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.



# copied from test_basic on Sep 13, 2022



# %matplotlib qt

# Go up by 2 directory and import 

import sys
import os.path as path
two_up =  path.abspath(path.join(__file__ ,"../.."))
sys.path.append(two_up)


# Actual imports 
import tmart
import numpy as np
from Py6S.Params.atmosprofile import AtmosProfile
import Py6S



# Specify wavelength in nm
band = Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_08)
band = None


wl = 400


### DEM and reflectance ###

# Three same-size numpy arrays are needed
image_DEM = np.array([[0,0],[0,0]]) # in meters
image_reflectance = np.array([[0.1,0.1],[0.1,0.1]]) # unitless     
image_isWater = np.zeros(image_DEM.shape) # 1 is water, 0 is land
image_isWater = np.array([[1,1],[1,1]])


# Synthesize a surface object
my_surface = tmart.Surface(DEM = image_DEM,
                           reflectance = image_reflectance,
                           isWater = image_isWater,
                           cell_size = 20_000)  

# Set background information, 1 or 2 background surfaces can be set;
# If 2: the first background is the one closer to [0,0]
my_surface.set_background(bg_ref        = [0.1,0.1], # background reflectance
                          bg_isWater    = [1,1], # if is water
                          bg_elevation  = 0, # elevation of both background
                          bg_coords     = [[0,0],[10,10]]) # a line dividing two background                                    


### Atmosphere ###

# Atmophere profile comes from 6S
atm_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeSummer) 
aerosol_type = 'Maritime' 
aot550 = 0

# Synthesize an atmosphere object     
my_atm = tmart.Atmosphere(atm_profile, aot550, aerosol_type)


### Running T-Mart ###

# Make a T-Mart object and modify wind and water parameters 
my_tmart = tmart.Tmart(Surface = my_surface, Atmosphere= my_atm, shadow=False)
my_tmart.set_wind(wind_speed=5, wind_dir=0)
my_tmart.set_water(water_salinity=35, water_temperature=20)


# Specify the sensor's position (x, y, z), viewing direction relative 
# to the sensor (zenith, azimuth), sun's direction relative to the target 
# (zenith, azimuth)
my_tmart.set_geometry(sensor_coords=[51,50,130_000], 
                      target_pt_direction=[180,0],
                      sun_dir=[0,0])

# Run



'''

results = my_tmart.run_plot(wl=wl, band=band, plot_on=True, plot_range=[0,100000,0,100000,0,100000])

'''


n_photon = 10_000
results = my_tmart.run(wl=wl, band=band, n_photon=n_photon,nc= 10,njobs= 100)
results = np.vstack(results)

# Calculate reflectances using recorded photon information 
R = tmart.calc_ref(results)
for k, v in R.items():
    print(k, '     ' , v)





''' Run the code above 


df = results 

# go to tm_calcref, run the assigning-column-name code 

# identify a photon with lots of movements in dfpd

df_output = pd.DataFrame()






'''
















