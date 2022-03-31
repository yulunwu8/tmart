# Other: DEM-shadow

Source: https://github.com/tomderuijter/python-dem-raycast

Modified by Yulun Wu, March 11, 2022

Slightly improved the computation effeciency by skiping the elevations 0s in the begining of the lines because they don't cast shadows, as often encountered in coastal areas. 



## Instruction 

First, import the libraries. Place the **python\_dem\_shadows** folder in the same folder as your script and data. 


```python
import os
import sys
import numpy as np
from osgeo import gdal, ogr, osr
import matplotlib.pyplot as plt
import python_dem_shadows as dem
print ("Imports done")
```

Modify your input:

```python
file_in = "cdem_16_17_26_19N.tif"
zenith = 55.54170463 
azimuth = 176.92407132
```


Load the DEM raster:

```python
ds = gdal.Open(file_in)

if ds is None:
    print ("Warning: can't open raster dataset")
    sys.exit() # error message

xsize = ds.RasterXSize
ysize = ds.RasterYSize
nbands = ds.RasterCount
projection = ds.GetProjection()
geotransform = ds.GetGeoTransform()

bands = []
for i in range(nbands):
    bands.append(ds.GetRasterBand(i+1))

arrays = []
for i in range(nbands):
    arrays.append(bands[i].ReadAsArray(0,0,xsize,ysize))

# Load DEM 
dem_data = arrays[0]
```

Plot the DEM in gray scale if needed:

```python
plt.imshow(dem_data, cmap='Greys_r')
plt.title("DEM")
plt.show()
```

Set up parameters:  

```python
# Solar angle vector 
sv = dem.normal_vector(zenith,  azimuth)

# cell_width 
dx = geotransform[1]

# cell height, abs because its value usually increases as it moves south 
dy = abs(geotransform[5])
```

Calculate shadow and plot it. In the output numpy array: 1 is no-shadow, 0 is shadow.


```python
# Cast shadow. Turn off print_progress if you don't need it
drop_shadow = dem.project_shadows(dem_data, sv, dx, dy, print_progress=True)

# Plot
plt.imshow(drop_shadow, cmap='Greys_r')
plt.title('Ray cast')
plt.show()
```

If needed, export the numpy array to a tiff file with the same geospatial information as the input raster.

```python
driver = gdal.GetDriverByName("GTiff")
dataset = driver.Create('test.tif', xsize, ysize, 1, gdal.GDT_Float32)
dataset.SetGeoTransform(geotransform)
dataset.SetProjection(projection)
band = dataset.GetRasterBand(1)
band.WriteArray(drop_shadow) 
dataset = None 
```












