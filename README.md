# T-Mart: Topography-adjusted Monte-carlo Adjacency-effect Radiative Transfer Code

## Description 

T-Mart solves the radiative transfer in a 3D surface-atmosphere system through a Monte Carlo approach. T-Mart features arbitrary surface models which allow simulating and correcting for the adjacency effect in aquatic remote sensing. 


## Links


Home page: <a href="https://github.com/yulunwu8/tmart" target="_blank">https://github.com/yulunwu8/tmart</a>

User guide: <a href="https://tmart-rtm.github.io" target="_blank">https://tmart-rtm.github.io</a>

## Publication

Wu, Y., Knudby, A., & Lapen, D. (2023). Topography-Adjusted Monte Carlo Simulation of the Adjacency Effect in Remote Sensing of Coastal and Inland Waters. *Journal of Quantitative Spectroscopy and Radiative Transfer*, 108589. <a href="https://doi.org/10.1016/j.jqsrt.2023.108589" target="_blank">https://doi.org/10.1016/j.jqsrt.2023.108589</a>

## Installation 

1 - Create a conda environment and activate it: 

```bash
conda create --name tmart python=3.9
conda activate tmart
```

2 - Install dependencies: 

```bash
conda install -c conda-forge Py6S
```


3 - Install tmart: 

```bash
pip3 install tmart
```

## Quick start: adjacency-effect correction 

T-Mart supports adjacency-effect correction for Sentinel-2 MSI and Landsat 8/9 OLI/OLI-2 products. Correction is performed directly on level-1 products and can be followed by any amtospheric correction tools. 

Minimal input: 

```python
import tmart
file = 'user/test/S2A_MSIL1C_20160812T143752_N0204_R096_T20MKB_20160812T143749.SAFE'

# NASA EarthData Credentials, OB.DAAC Data Access needs to be approved
username = 'abcdef'
password = '123456'

### Multiprocessing needs to be wrapped in 'if __name__ == "__main__":' for Windows systems, this is optional for Mac OS
if __name__ == "__main__":
    tmart.AEC.run(file, username, password)
```

See <a href="https://tmart-rtm.github.io/ins_aec.html" target="_blank">Instruction - Adjacency-Effect Correction</a> for detailed instructions.

## Quick start: adjacency-effect modelling

```python
import tmart
import numpy as np
from Py6S.Params.atmosprofile import AtmosProfile

# Specify wavelength in nm
wl = 400

### DEM and reflectance ###
image_DEM = np.array([[0,0],[0,0]]) # in meters
image_reflectance = np.array([[0.1,0.1],[0.1,0.1]]) # unitless     
image_isWater = np.zeros(image_DEM.shape) # 1 is water, 0 is land

# Synthesize a surface object
my_surface = tmart.Surface(DEM = image_DEM,
                           reflectance = image_reflectance,
                           isWater = image_isWater,
                           cell_size = 10_000)  
                               
### Atmosphere ###
atm_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeSummer) 
my_atm = tmart.Atmosphere(atm_profile, aot550 = 0, aerosol_type = 'Maritime')

### Create T-Mart Object ###
my_tmart = tmart.Tmart(Surface = my_surface, Atmosphere= my_atm, shadow=False)
my_tmart.set_geometry(sensor_coords=[51,50,130_000], # x, y, and z 
                      target_pt_direction=[180,0], # zenith and azimuth
                      sun_dir=[0,0]) # zenith and azimuth

### Multiprocessing needs to be wrapped in 'if __name__ == "__main__":' for Windows systems. 
### This can be skipped for Unix-based systems. 
if __name__ == "__main__":
    results = my_tmart.run(wl=wl, band=None, n_photon=10_000)
    
    # Calculate reflectances using recorded photon information 
    R = tmart.calc_ref(results)
    for k, v in R.items():
        print(k, '     ' , v)

```

The output should be similar to this: 

```
========= Initiating T-Mart =========
Number of photons: 10000
Using 10 core(s)
Number of job(s): 100
Wavelength: 400
target_pt_direction: [180, 0]
sun_dir: [0, 0]
=====================================
Jobs remaining = 102
Jobs remaining = 72
Jobs remaining = 42
Jobs remaining = 12
=====================================
Calculating radiative properties...
R_atm       0.12760589889823587
R_dir       0.06046419017201067
R_env       0.012888590547129805
R_total       0.20095867961737635

```

See user guide for more detailed instructions. 


## Known issue(s)

`rasterio` version 1.3.11 leads to unprojected S2 image files when used in Mac’s Terminal. 

- Workarounds: run the code in an IDE, or downgrade to `rasterio` version 1.3.9 for terminal-based workflows.


## Others

T-Mart can calculate reflectances of various units, see Table 1 in Wu et al. (2023) for examples. 

For questions and suggestions (which I'm always open to!), please email Yulun at [yulunwu8@gmail.com](mailto:yulunwu8@gmail.com)

