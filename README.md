# T-Mart: Topography-adjusted Monte-carlo Adjacency-effect Radiative Transfer Code

## Description 

T-Mart solves the radiative transfer in a 3D ocean-atmosphere system through a Monte-Carlo approach. T-Mart features arbitrary surface models which allow simulations of the adjacnecy effect in aquatic remote sensing. 


## Links


Home page: <a href="https://github.com/yulunwu8/tmart" target="_blank">https://github.com/yulunwu8/tmart</a>

User guide: <a href="https://tmart-rtm.github.io" target="_blank">https://tmart-rtm.github.io</a>

## Installation 

1 - Create a conda environment and activate it: 

```
conda create --name tmart python=3.9
conda activate tmart
```

2 - Install dependencies: 

```
conda install -c conda-forge Py6S numpy pandas scipy pathos matplotlib
```

3 - Install tmart: 

```
pip3 install tmart
```

## Test Run

```
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
results = np.vstack(results)

# Calculate reflectances using recorded photon information 
R = tmart.calc_ref(results)
for k, v in R.items():
    print(k, '     ' , v)

```

Output should be similar to below: 

```
========= Initiating T-Mart =========
Number of photons: 10000
Using 10 core(s)
Number of job(s): 100
Wavelength: 400
target_pt_direction: [180, 0]
sun_dir: [0, 0]
=====================================
Tasks remaining = 102
Tasks remaining = 72
Tasks remaining = 42
Tasks remaining = 12
=====================================
Calculating reflectances...
R_atm       0.12576169910706145
R_dir       0.11418041191253057
R_env       0.013670152459801528
R_total       0.25361226347939353

```



## Other


For questions and suggestions (which I'm always open to!), please email [yulunwu8@gmail.com](mailto:yulunwu8@gmail.com)