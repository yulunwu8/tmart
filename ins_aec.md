# Instruction - Adjacency-Effect Correction

## Minimal Input

A single function is used to perform adjacency-effect correction (AEC) in T-Mart. Correction is performed directly on level-1 products, and the output is adjacency-effect-free top-of-atmosphere products in the same format as the input level-1 products, therefore this workflow can be followed by any amtospheric-correction tools. Currently it only supports Sentinel-2 MSI and Landsat 8 OLI imagery. 

NASA EarthData Credentials are needed to retrieve ozone, water vapour, and aerosol information for accurate AEC. You may need to approve OB.DAAC Data Access in your <a href="https://urs.earthdata.nasa.gov/profile" target="_blank">EarthData account</a>.


Minimal input to the AEC.run function includes path to satellite files and EarthData Credentials. See the <a href="https://tmart-rtm.github.io/tmart.html#module-tmart.AEC.run" target="_blank">AEC.run Function</a> tab for all arguments. 

```python
file = 'user/test/S2A_MSIL1C_20160812T143752_N0204_R096_T20MKB_20160812T143749.SAFE'
username = 'abcdef'
password = '123456'

### Multiprocessing needs to be wrapped in 'if __name__ == "__main__":' for Windows systems
if __name__ == "__main__":
    tmart.AEC.run(file, username, password)
```

## Overwrite Existing Files

By default, this creates a copy of the original satellite files in the same directory, in a new folder that starts with 'AEC_'. To overwrite the eixsting files, add 'overwrite=True' to the arguments: 

```python
tmart.AEC.run(file, username, password, overwrite=True)
```

## Ancillary and Log Files  

During the AEC process, a number of files are generated: 

- **tmart\_log\_\*.txt**: detailed processing information, as printed in the Python console. 
- **tmart\_atm\_info.txt**: atmosphere and aerosol information used in the processing. This includes aerosol type, angstrom exponent, single scattering albedo, AOT at 550 nm, total column ozone, and total precipitable water vapour. 
- **tmart\_ancillary/\*.nc**: ancillary files from NASA Ocean Color. 
- **tmart\_completed.txt**: a record of bands that have been corrected for the adjacency effect.  

## AEC Configuration

A configuration file is stored in the *tmart* package folder. Brief descriptions are given in the file. Most of the configuration settings are tuned for best performance. In case a large amount of water pixels are falsely masked as land, ``AE_land`` can be set as True in order to perform AEC across the entire scene. 

## Additional Arguments 

Lastly, ``AOT`` and ``n_photon`` can be specified. You can specify ``AOT`` if you are certain about its value, and this reduces processing time. ``n_photon`` is the number of photons used in each T-Mart run, 100_000 is recommended for accurate results. It can be reduced to 10_000 for quicker computation. 

```python
tmart.AEC.run(file, username, password, overwrite=True, AOT = 0.05, n_photon = 10_000)
```

## Assumptions in the Processing 

A few assumptions are made in the processing, violations can lead to various degrees of errors in the output AE-free product: 

- Isotropic/Lambertian surface 
- Homogeneous atmosphere and aerosols across the scene 
- Flat surface or lack of topography
- Reflectance outside the scene is the median of the scene at each wavelength, it is homogeneous and it extends to infinity









