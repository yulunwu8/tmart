# Overview

## T-Mart: Topography-adjusted Monte-carlo Adjacency-effect Radiative Transfer Code


This code models the radiative transfer in a 3D atmosphere-ocean system. In addition to the Monte Carlo-based radiative transfer functions, there are three environmental components in the code: atmosphere, water and land. 

- The atmosphere consists of layers with various scattering and absorbing properties.  
- Land is assumed to be Lambertian. The topography of land is modelled by triangulating the pixels of an input DEM. 
- Water has three reflectance properties: 1) water-leaving reflectance, 2) white-caps and 3) surface reflectance (sky glint and sun glint). 1 and 2 are assumed to be Lambertian, 3 is calculated through Cox&Munk slope statistics. 

Yulun Wu | March 13, 2022 | Requested by Liquid Geomatics | For any questions, please email ywu146@uottawa.ca

## Required Libraries

T-Mart is written in Python 3, and it requires: Py6S, matplotlib, pathos, multiprocessing, numpy, pandas.


## Input

**Essential**

- Wavelength: wavelength in nm.
- Surface: 
	- DEM: Digital Elevation Model, the elevation of pixels in this run.
	- Cell size: the width and length of each pixel. 
	- Reflectance: reflectance of land or water-leaving reflectance of water, Lambertian. 
	- is_water: which pixels are water pixels. 
	- Background information: the elevation and reflectance of the background beyond the range of the pixels.

- Atmosphere profile: choose from Py6S. 

- n_photon: number of photons used in each run, default 1000. Recommend 10,000 to 100,000.
- Geometry: photon starting position, solar angle, viewing angle. 

**Geometry**
 
 The geometry in T-Mart follows the diagram below.
 
 
![Geometry](files/geometry.png)


**Optional**

- Atmosphere 
 	- AOT550: aerosol optical thickness at 550nm.
	- Aerosol_SPF: normalized aerosol scattering phase function.
	- Aerosol scale height: default 2 km.
	- n_layers: number of atmospheric layers: default 10. Having more layers may slightly increase the computation time
- Wind speed
- Water salinity: in unit of parts per thousand, default 0.
- Water temperature: in Celsius, default 25.




## Output 

Reflectances (definitions following 6S):

- Atmospheric intrinsic reflectance
- Environmental reflectance 
- Direct reflectance









