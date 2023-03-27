


import tmart
import numpy as np
from Py6S.Params.atmosprofile import AtmosProfile

# Specify wavelength in nm
wl = 400

# DEM and reflectance ###
image_DEM = np.array([[0,0],[0,0]]) # in meters
image_reflectance = np.array([[0.1,0.1],[0.1,0.1]]) # unitless     
image_isWater = np.zeros(image_DEM.shape) # 1 is water, 0 is land
image_isWater = np.array([[1,1],[1,1]])


# Synthesize a surface object
my_surface = tmart.Surface(DEM = image_DEM,
                           reflectance = image_reflectance,
                           isWater = image_isWater,
                           cell_size = 10_000)  
                               
### Atmosphere ###
atm_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeSummer) 
aerosol_type = 'Maritime'  
my_atm = tmart.Atmosphere(atm_profile, aot550 = 0, aerosol_type = 'Maritime'  )

### Running T-Mart ###
my_tmart = tmart.Tmart(Surface = my_surface, Atmosphere= my_atm, shadow=False)
my_tmart.set_geometry(sensor_coords=[51,50,130_000], 
                      target_pt_direction=[180,0],
                      sun_dir=[0,0])

results = my_tmart.run(wl=wl, band=None, n_photon=10_000,nc= 10,njobs= 100)
# results = np.vstack(results)

# Calculate reflectances using recorded photon information 
R = tmart.calc_ref(results)
for k, v in R.items():
    print(k, '     ' , v)
    
    
'''
    

import tmart
import numpy as np
from Py6S.Params.atmosprofile import AtmosProfile

# Specify wavelength in nm
wl = 800

### DEM and reflectance ###

# Three same-size numpy arrays are needed
image_DEM = np.array([[30_000,20_000],[5_000,10_000]]) # in meters
image_reflectance = np.array([[0.5,0.5],[0.5,0.5]])  # unitless     
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
my_atm = tmart.Atmosphere(atm_profile)

### Running T-Mart ###

# Make a T-Mart object 
my_tmart = tmart.Tmart(Surface = my_surface, Atmosphere= my_atm)

# Specify the sensor's position (x, y, z), viewing direction relative 
# to the sensor (zenith, azimuth), sun's direction relative to the target 
# (zenith, azimuth)
my_tmart.set_geometry(sensor_coords=[18_000,17_000,130_000], 
                      target_pt_direction=[180,0],
                      sun_dir=[0,0])

# Run and plot on a separate window
%matplotlib qt
results = my_tmart.run_plot(wl=wl, plot_on=True, plot_range=[0,100_000,0,100_000,0,100_000])
results = np.vstack(results)

# Calculate reflectances using recorded photon information 
R = tmart.calc_ref(results)
for k, v in R.items():
    print(k, '     ' , v)


'''



