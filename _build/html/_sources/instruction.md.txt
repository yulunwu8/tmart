# Instruction

We start with a simple setting. I will introduce new functions as we go, and have a 'put together' code in each sub-section.

To run the code below, you can overwrite the example scripts in the main folder, create your own script in the main folder, or copy the *./tmart* folder to where your script is. 


## Quick Start

Run T-Mart with simple settings here. 


```python
import tmart
import numpy as np
from Py6S.Params.atmosprofile import AtmosProfile

# Specify wavelength in nm
wl = 400


### DEM and reflectance ###

# Three same-size numpy arrays are needed
image_DEM = np.array([[0,0],[0,0]]) # in meters
image_reflectance = np.array([[0.1,0.1],[0.1,0.1]]) # unitless     
image_isWater = np.zeros(image_DEM.shape) # 1 is water, 0 is land

# pixel width and length, in meters 
cell_size = 20_000 

# Synthesize a surface object
my_surface = tmart.Surface(image_DEM,image_reflectance,image_isWater,cell_size)  


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
my_tmart.set_geometry(sensor_coords=[51,50,130_000], 
                      target_pt_direction=[180,0],
                      sun_dir=[0,0])

results = my_tmart.run(wl=wl)
results = np.vstack(results)

# Calculate reflectances using recorded photon information 
R = tmart.calc_ref(results)
print(R)
```


## Multiple Processing 

Monte Carlo simulations have inherent noise - it decreases with a larger sample size. Multiprocessing is used to speed up the computation. 

By default, 10,000 photons are used. Number of CPU cores to use is automated. Default number of jobs is 80. Use 100,000 photons for more stable results (it may take several minutes for a single computation).



```python
# Number of photons
n_photon = 100_000

# Number of CPU cores to use, sometimes default doesn't work and 
# you need to specify the number of cores you have. 
nc = 10

# Dividing the task into n jobs, I found 80-120 function 
# more or less the same. 
njobs = 100

results = my_tmart.run(wl=wl, n_photon=n_photon,nc= nc,njobs= njobs)
```

Now we add the snippet to the code and put things together:

```python
import tmart
import numpy as np
from Py6S.Params.atmosprofile import AtmosProfile

# Specify wavelength in nm
wl = 400


### DEM and reflectance ###

# Three same-size numpy arrays are needed
image_DEM = np.array([[0,0],[0,0]]) # in meters
image_reflectance = np.array([[0.1,0.1],[0.1,0.1]]) # unitless     
image_isWater = np.zeros(image_DEM.shape) # 1 is water, 0 is land

# pixel width and length, in meters 
cell_size = 20_000 

# Synthesize a surface object
my_surface = tmart.Surface(image_DEM,image_reflectance,image_isWater,cell_size)  


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
my_tmart.set_geometry(sensor_coords=[51,50,130_000], 
                      target_pt_direction=[180,0],
                      sun_dir=[0,0])

# Run
n_photon = 10_000
nc = 10
njobs = 100
results = my_tmart.run(wl=wl, n_photon=n_photon,nc= nc,njobs= njobs)
results = np.vstack(results)

# Calculate reflectances using recorded photon information 
R = tmart.calc_ref(results)
print(R)
```

## Observing the Movements of a Single Photon

Instead of running lots of photons, we can run a single photon and observe where it goes and what happens. This is mostly for debugging purposes. The details of each movement will be printed. 

We replace Tmart's *run* function with the *run_plot* function here, note that number of photons and multiprocessing are not needed here. 

```python
results = my_tmart.run_plot(wl=wl)
```

We can also plot the photon's movements by turnining on *plot_on*: 

```python
results = my_tmart.run_plot(wl=wl, plot_on=True)
```


By default, the limits of X, Y, Z are all 0 to 100,000. We can specify them in the format of [xmin, xmax, ymin, ymax, zmin, zmax]:

```python
results = my_tmart.run_plot(wl=wl, plot_on=True, plot_range=[0,100_000,0,100_000,0,100_000])
```

Lastly, to switch to an interactive plotting mode where you can zoom and rotate the camera in a plot interactively, and switch back, use: 

```python
# Interactive mode
%matplotlib qt

# Inline plotting 
%matplotlib inline
```

Putting it all together. We can still calculate reflectances using one photon thanks to the local-estimate technique, but this is not accurate - precise stats come from a large sample size. 

```python
import tmart
import numpy as np
from Py6S.Params.atmosprofile import AtmosProfile


# Specify wavelength in nm
wl = 400


### DEM and reflectance ###

# Three same-size numpy arrays are needed
image_DEM = np.array([[0,0],[0,0]]) # in meters
image_reflectance = np.array([[0.1,0.1],[0.1,0.1]]) # unitless     
image_isWater = np.zeros(image_DEM.shape) # 1 is water, 0 is land

# pixel width and length, in meters 
cell_size = 20_000 

# Synthesize a surface object
my_surface = tmart.Surface(image_DEM,image_reflectance,image_isWater,cell_size)  


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
my_tmart.set_geometry(sensor_coords=[51,50,130_000], 
                      target_pt_direction=[180,0],
                      sun_dir=[0,0])

# Run and plot on a separate window
%matplotlib qt
results = my_tmart.run_plot(wl=wl, plot_on=True, plot_range=[0,100_000,0,100_000,0,100_000])
results = np.vstack(results)

# Calculate reflectances using recorded photon information 
R = tmart.calc_ref(results)
print(R)
```
Here's the first plot generated by the code above, the following ones are not shown here.
![Geometry](files/movement.png)

There are lines of different colours in the plot: 

- Blue to red: the photon moves from blue to red.
- Green: surface normal of a Lambertian surface.
- Dark blue: surface normal of a specular surface.
- Orange: next moving direction 

## Water Pixels 

From here we go back to multiprocessing lots of photons. 

We can change the pixels from land to water by specifying the water pixels in the image_isWater numpy array:

```python
image_isWater = np.array([[1,1],[1,1]])
```




## Modify Background 

By default, the background surface takes the average reflectance and average elevation of the DEM and the reflecting surface. 

A total of two background surfaces can be specified to model coastal environments, they are divided by a line connected by two specified coordinates (*bg_coords*). We can modify the reflectance and if-is-water of each of the two background surfaces. We can modify the elevation of the background too, but both background surfaces have to share the same elevation (to avoid gaps in between). 


```python
# Set background information, 1 or 2 background surfaces can be set;
# If 2: the first background is the one closer to [0,0]
my_surface.set_background(bg_ref        = [0.02,0.02], # background reflectance
                          bg_isWater    = [0,0], # if is water
                          bg_elevation  = 0, # elevation of both background
                          bg_coords     = [[0,0],[10,10]]) # a line dividing two background                                    
```

## Adding Aerosol 

Add AOT550 and aerosol scattering phase function into the atmosphere object. Currently T-Mart only supports the maritime aerosol model from 6S. If you need to use other models, let me know. The miepython package is handy — I used 6S' input and got the same scattering properties as 6S thanks to help from John Hedley. I just haven't had the incentive to look into it...

```python
# Atmosphere
atm_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeSummer) 
aerosol_SPF = 'tmart/ancillary/aerosol_maritime_SPF.csv' 
aot550 = 0.1
my_atm = tmart.Atmosphere(atm_profile, aot550, aerosol_SPF)
```



## Additional Parameters 

Set the number of atmosphere layers and aerosol scale height. T-Mart can calculate the movement of a photon through multiple layers so having lots of layers (e.g., 20) does not slow down the computation significantly. 

```python
n_layers = 20
aerosol_scale_height = 2 # Unless you have a reason, don't change this
    
my_atm = tmart.Atmosphere(atm_profile, aot550, aerosol_SPF, n_layers, aerosol_scale_height)
```
Once we have the Tmart object, we can modify wind and water properties. They are both related to the specular reflectance of the water surface. 

```python
my_tmart = tmart.Tmart(Surface = my_surface, Atmosphere= my_atm)
my_tmart.set_wind(wind_speed=5, wind_dir=0)
my_tmart.set_water(water_salinity=35, water_temperature=20)
```


## Full Single Run 

Now we can put everything together. Below is the code I usually use for my single runs (different components of TOA reflectance at a single wavelength in a single scenario, but with lots of photons):

```python
import tmart
import numpy as np
from Py6S.Params.atmosprofile import AtmosProfile

# Specify wavelength in nm
wl = 400


### DEM and reflectance ###

# Three same-size numpy arrays are needed
image_DEM = np.array([[0,0],[0,0]]) # in meters
image_reflectance = np.array([[0.1,0.1],[0.1,0.1]]) # unitless     
image_isWater = np.zeros(image_DEM.shape) # 1 is water, 0 is land

# pixel width and length, in meters 
cell_size = 20_000 

# Synthesize a surface object
my_surface = tmart.Surface(image_DEM,image_reflectance,image_isWater,cell_size)  
# Set background information, 1 or 2 background surfaces can be set;
# If 2: the first background is the one closer to [0,0]
my_surface.set_background(bg_ref        = [0.1,0.1], # background reflectance
                          bg_isWater    = [1,1], # if is water
                          bg_elevation  = 0, # elevation of both background
                          bg_coords     = [[0,0],[10,10]]) # a line dividing two background                                    

### Atmosphere ###

# Atmophere profile comes from 6S
atm_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeSummer) 
aerosol_SPF = 'tmart/ancillary/aerosol_maritime_SPF.csv' 
aot550 = 0.1
n_layers = 20
aerosol_scale_height = 2 # Unless you have a reason, don't change this

# Synthesize an atmosphere object     
my_atm = tmart.Atmosphere(atm_profile, aot550, aerosol_SPF, n_layers, aerosol_scale_height)


### Running T-Mart ###

# Make a T-Mart object 
my_tmart = tmart.Tmart(Surface = my_surface, Atmosphere= my_atm)
my_tmart.set_wind(wind_speed=5, wind_dir=0)
my_tmart.set_water(water_salinity=35, water_temperature=20)

# Specify the sensor's position (x, y, z), viewing direction relative 
# to the sensor (zenith, azimuth), sun's direction relative to the target 
# (zenith, azimuth)
my_tmart.set_geometry(sensor_coords=[51,50,130_000], 
                      target_pt_direction=[180,0],
                      sun_dir=[0,0])

# Run
n_photon = 10_000
nc = 10
njobs = 100
results = my_tmart.run(wl=wl, n_photon=n_photon,nc= nc,njobs= njobs)
results = np.vstack(results)

# Calculate reflectances using recorded photon information 
R = tmart.calc_ref(results)
print(R)
```


## Be Creative at What You Can Do! 

You can do lots of things with T-Mart! E.g.:

- Load your own DEM and image (they have to be the same size),
- Loop through wavelengths, distances from shore, AOT, viewing or solar angles, etc. 

T-Mart currently provides four built-in spectral landcover types: soil, vegetation, water and water with a chlorophyll concentration of 1 µg/L. For example: 


```python
# Create object
water = tmart.spectral_surface('water_chl1')
      
# Find spectral reflectance at 400nm
water.wl(400)
```

## Case Study: Typical Ocean-Colour Observation 

Here we model a typical hyperspectral ocean-colour observation with the following parameters: 

- View zenith = 30°
- Solar zenith = 30°
- View azimuthal = 90° (to minimize sun-glint)
- Windspeed = 5 m/s
- Atmosphere profile: mid-latitude summer
- Aerosol profile: maritime aerosols 
- AOT550: 0.1 
- Water salinity: 35 parts per thousand
- Water temperature: 25 C
- Water chl-a concentration: 1 µg/L

Water's specular reflectance is calculated following Cox&Munk slope statistics. 

We loop through wavelengths 400 to 1100nm at an interval of 10nm, use a Python list to collect the reflectances at each wavelength, convert it to a Pandas dataframe, and plot the reflectances with the *matplotlib* package:

```python
import tmart
import numpy as np
import pandas as pd
from Py6S.Params.atmosprofile import AtmosProfile

water = tmart.spectral_surface('water_chl1')

# Specify wavelength in nm

pd_output = []

wavelengths = range(400,1110,10)

for wl in wavelengths:

    ### DEM and reflectance ###
    
    # Three same-size numpy arrays are needed
    image_DEM = np.array([[0,0],[0,0]]) # in meters
    image_reflectance = np.array([[water.wl(wl),water.wl(wl)],
                                  [water.wl(wl),water.wl(wl)]]) # unitless     
    image_isWater = np.array([[1,1],[1,1]]) # 1 is water, 0 is land
    
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
    aerosol_SPF = 'tmart/ancillary/aerosol_maritime_SPF.csv' 
    aot550 = 0.1
    n_layers = 20
    aerosol_scale_height = 2 # Unless you have a reason, don't change this
    
    # Synthesize an atmosphere object     
    my_atm = tmart.Atmosphere(atm_profile, aot550, aerosol_SPF, n_layers, aerosol_scale_height)
    
    
    ### Running T-Mart ###
    
    # Make a T-Mart object 
    my_tmart = tmart.Tmart(Surface = my_surface, Atmosphere= my_atm)
    my_tmart.set_wind(wind_speed=5, wind_dir=0)
    my_tmart.set_water(water_salinity=35, water_temperature=20)
    
    # Specify the sensor's position (x, y, z), viewing direction relative 
    # to the sensor (zenith, azimuth), sun's direction relative to the target 
    # (zenith, azimuth)
    my_tmart.set_geometry(sensor_coords=[51,50,130_000], 
                          target_pt_direction=[150,90],
                          sun_dir=[30,0])
    
    # Run
    n_photon = 10_000
    nc = 10
    njobs = 100
    results = my_tmart.run(wl=wl, n_photon=n_photon,nc= nc,njobs= njobs)
    results = np.vstack(results)
    
    # Calculate reflectances using recorded photon information 
    R = tmart.calc_ref(results)
    print(R)
    R['Wavelength'] = wl
    
    pd_output.append(R)

pd_output = pd.DataFrame(pd_output)

# Plot 
import matplotlib.pyplot as plt

plt.plot(pd_output["Wavelength"], pd_output["R_atm"], color='grey', label='Atmospheric intrinsic R.')
plt.plot(pd_output["Wavelength"], pd_output["R_dir"], color='blue', label='Direct R.')
plt.plot(pd_output["Wavelength"], pd_output["R_env"], color='orange', label='Environmental R.')
plt.plot(pd_output["Wavelength"], pd_output["R_total"], color='black', label='Total R.')

plt.legend()
plt.show()
```
![waterTOA](files/waterTOA.png)

## Useful Information 

The following two pandas dataframes completely describe the atmosphere in T-Mart. 



```python
# The atmosphere profile
my_atm_profile = my_tmart.atm_profile_wl

# The aerosol scattering phase function
my_aerosol_SPF = my_tmart.aerosol_SPF_wl
```

The columns of the atmosphere profile:

- Alt_bottom: the bottom the the layer in km.
- Alt_top: the top of the layer in km.
- ot_abs: optical thickness of absroption. 
- ot_rayleigh: optical thickness of molecular or rayleigh scattering. No scattering phase function needed because it is built-in.
- ot_mie: optical thickness of aerosol scattering, relies on the input scattering phase function (see below).
- ot_scatt: sum of ot_rayleigh and ot_mie. 

The columns of the aerosol scattering phase function (currently rely on pre-calculated values, may use miepython in the future to compute on the go):

- Angle: the scattering angle. 
- Value: normalized intensity at the corresponding scattering angle. 








