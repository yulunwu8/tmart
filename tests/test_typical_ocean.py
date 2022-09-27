import sys
import os.path as path
two_up =  path.abspath(path.join(__file__ ,"../.."))
sys.path.append(two_up)






import tmart
import numpy as np
import pandas as pd
from Py6S.Params.atmosprofile import AtmosProfile

water = tmart.SpectralSurface('water_chl1')

# Master output
pd_output = []

# Specify wavelength in nm
wavelengths = range(400,1110,10)

for wl in wavelengths:

    ### DEM and reflectance ###
    
    # Three same-size numpy arrays are needed
    image_DEM = np.full((2, 2), 0) # in meters
    image_reflectance = np.full((2, 2), water.wl(wl)) # unitless       
    image_isWater = np.full((2, 2), 0) # 1 is water, 0 is land
    
    # pixel width and length, in meters 
    cell_size = 20_000 
    
    # Synthesize a surface object
    my_surface = tmart.Surface(image_DEM,image_reflectance,image_isWater,cell_size)  
    # Set background information, 1 or 2 background surfaces can be set;
    # If 2: the first background is the one closer to [0,0]
    my_surface.set_background(bg_ref        = [water.wl(wl),water.wl(wl)], 
                              bg_isWater    = [1,1], # water
                              bg_elevation  = 0, # elevation of both background
                              bg_coords     = [[0,0],[10,10]]) # a line dividing two background                                    
    
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
    my_tmart = tmart.Tmart(Surface = my_surface, Atmosphere= my_atm, shadow=False)
    my_tmart.set_wind(wind_speed=10, wind_azi_avg = True)
    my_tmart.set_water(water_salinity=35, water_temperature=20)
    
    # Specify the sensor's position (x, y, z), viewing direction relative 
    # to the sensor (zenith, azimuth), sun's direction relative to the target 
    # (zenith, azimuth)
    my_tmart.set_geometry(sensor_coords=[51,50,130_000], 
                          target_pt_direction=[150,90],
                          sun_dir=[30,0])
    
    # Run
    n_photon = 10_000
    results = my_tmart.run(wl=wl, n_photon=n_photon)
    
    # Calculate reflectances using recorded photon information 
    R = tmart.calc_ref(results)
    
    for k, v in R.items():
        print(k, '     ' , v)
        
    # Add wavelength and append to output 
    R['Wavelength'] = wl
    pd_output.append(R)

pd_output = pd.DataFrame(pd_output)

# Plot 
import matplotlib.pyplot as plt

plt.plot(pd_output["Wavelength"], pd_output["R_atm"], color='grey', label='Atmospheric intrinsic R.')
plt.plot(pd_output["Wavelength"], pd_output["R_dir"], color='blue', label='Direct R.')
plt.plot(pd_output["Wavelength"], pd_output["R_env"], color='orange', label='Environmental R.')
plt.plot(pd_output["Wavelength"], pd_output["R_total"], color='black', label='Total R.')
plt.ylim(ymin=0)

plt.legend()
plt.show()




