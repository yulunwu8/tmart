# This file is part of TMart.
#
# Copyright 2024 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# Extract ozone and water vapour from NASA Ocean Color ancillary data 
def anci_get_OWV(files, metadata): 
    
    print('\nRetrieving ozone and water vapour at center of image: ')
    OWV_0 = _get_OWV(files[0], metadata['lat'], metadata['lon'])
    OWV_1 = _get_OWV(files[1], metadata['lat'], metadata['lon'])
    
    # Interpolate aerosol linearly in time
    minute = float(metadata['time'][14:16])
    second = float(metadata['time'][17:19])
    ratio = minute/60 + second/3600
    
    OWV = OWV_1 * ratio + OWV_0 * (1-ratio)
    print('Interpolated total column ozone is {:.2f} DU, total precipitable water vapour is {:.2f} kg m-2'.format(OWV[0],OWV[1]))
    return OWV
    
# Find value of ancillary data at lat and lon
def _get_OWV(file, lat, lon):
    
    import netCDF4 as nc4
    import numpy as np
    
    from scipy import interpolate
    
    # from https://stackoverflow.com/questions/31820107/is-there-a-numpy-function-that-allows-you-to-specify-start-step-and-number
    def linspace(start, stop, step=1.):
        import numpy as np
        """
        Like np.linspace but uses step instead of num
        This is inclusive to stop, so if start=1, stop=3, step=0.5
        Output is: array([1., 1.5, 2., 2.5, 3.])
        """
        return np.linspace(start, stop, int((stop - start) / step + 1))
    
    nc = nc4.Dataset(file, 'r') 
    
    lons = linspace(nc.getncattr('geospatial_lon_min'),
                    nc.getncattr('geospatial_lon_max'),
                    nc.getncattr('geospatial_lon_resolution'))
    lats = linspace(nc.getncattr('geospatial_lat_min'),
                    nc.getncattr('geospatial_lat_max'),
                    nc.getncattr('geospatial_lat_resolution'))
    
    ozone = nc['TO3'][:]
    WV = nc['TQV'][:]
    
    # Ozone
    interp_ozone = interpolate.RegularGridInterpolator((lats,lons), ozone)
    i_ozone = interp_ozone([lat,lon])[0]
    
    # Water vapour
    interp_WV = interpolate.RegularGridInterpolator((lats,lons), WV)
    i_WV = interp_WV([lat,lon])[0]
    
    return np.array([i_ozone, i_WV])

