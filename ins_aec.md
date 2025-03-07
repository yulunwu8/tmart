# Instruction - Adjacency-Effect Correction

## Mininum input

The *AEC.run* function is used to perform adjacency-effect correction (AEC) in T-Mart. Correction is performed directly on level-1 products, and the output is adjacency-effect-free top-of-atmosphere products in the same format as the input level-1 products, therefore this workflow can be followed by any amtospheric correction tools. Currently it supports Sentinel-2 MSI, Landsat 8/9 OLI/OLI-2 and PRISMA imagery.

NASA EarthData credentials are needed to retrieve ozone, water vapour, and aerosol optical thickness and composition for accurate AEC. You may need to approve OB.DAAC Data Access in your <a href="https://urs.earthdata.nasa.gov/profile" target="_blank">EarthData account</a>: click *Authorized Apps* under the *Applications* tab - click the *APPROVE MORE APPLICATIONS* button near the bottom, and authorize *OB.DAAC Data Access*.


Minimum input to the *AEC.run* function includes a path to satellite files and EarthData credentials. See <a href="https://tmart-rtm.github.io/tmart.html#module-tmart.AEC.run" target="_blank">AEC.run Function</a> for all the arguments. 

```python
import tmart
file = 'user/test/S2A_MSIL1C_20160812T143752_N0204_R096_T20MKB_20160812T143749.SAFE'
username = 'abcdef'
password = '123456'

# T-Mart uses multiprocessing, which needs to be wrapped in 'if __name__ == "__main__":' for Windows users. This is optional for Unix-based systems
if __name__ == "__main__":
    tmart.AEC.run(file, username, password)
```

## Overwrite existing files

By default, the processing creates a copy of the original satellite files in the same directory, in a new folder named 'AEC_*'. To overwrite the eixsting files, add 'overwrite=True' to the arguments: 

```python
tmart.AEC.run(file, username, password, overwrite=True)
```

## Ancillary and log files  

During the AEC process, a number of files are generated: 

- **tmart\_log\_\*.txt**: detailed processing information, as printed in the Python console. 
- **tmart\_atm\_info\_\*.txt**: atmosphere and aerosol information used in the processing. This includes aerosol type, angstrom exponent, single scattering albedo, AOT at 550 nm, total column ozone, and total precipitable water vapour. 
- **tmart\_ancillary/\*.nc**: GMAO MERRA2 ancillary files from the NASA Ocean Biology Processing Group. 
- **tmart\_completed.txt**: a record of bands that have been corrected for the adjacency effect. 
- **tmart\_preview.txt**: a image preview of pixels identified as water (only run when ``AE_land`` is False).

## AEC configuration

A TXT configuration file is stored in the *tmart* package folder. Its path is printed in the Python console and the log file. Brief descriptions are given in the file. Most of the configuration settings are tuned for best performance. 

By default, T-Mart identifies water pixels and only modify their values, leaving land pixel values unchanged to facilitate the existing calibration of atmospheric correction processors that extract information from land pixels. In case a significant number of water pixels are falsely masked as land, the ``mask_SWIR_threshold`` (the reflectance threshold in a SWIR band used to mask non-water pixels; default value: 0.03) can be increased based on the water pixel values in the scene. Modifying ``mask_SWIR_threshold`` in the *AEC.run* function overwrites the value in config.txt. Alternatively, setting ``AE_land`` to True enables AEC across the entire scene. 

## Additional arguments 

``AOT`` and ``n_photon`` can be specified manually. ``AOT`` is the aerosol optical thickness at 550 nm, you can specify it if you are certain about its value or simply to test the impact of using different values. ``n_photon`` is the number of photons used in each T-Mart run; the default value of 100,000 is recommended for accurate results. It can be reduced to 10,000 for quicker computation. 

```python
tmart.AEC.run(file, username, password, overwrite=True, AOT = 0.05, n_photon = 10_000)
```

## Assumptions in AEC 

A few assumptions are made in the processing, violations can lead to various degrees of errors in the output AE-free product: 

- Isotropic/Lambertian surface 
- Vertically stratified but horizontally homogeneous atmospheric molecules and aerosols across the scene 
- Flat surface or lack of topography

See <a href="https://tmart-rtm.github.io/tmart.html#module-tmart.AEC.run" target="_blank">AEC.run Function</a> for all the arguments. 







