# This file is part of TMart.
#
# Copyright 2024 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# Extract aerosol type from NASA Ocean Color ancillary data 
def anci_get_AER(files, metadata): 
    
    print('\nRetrieving aerosol type at center of image: ')
    r_maritime_0 = _get_AER(files[0], metadata['lat'], metadata['lon'])
    r_maritime_1 = _get_AER(files[1], metadata['lat'], metadata['lon'])
    
    # Interpolate aerosol linearly in time
    minute = float(metadata['time'][14:16])
    second = float(metadata['time'][17:19])
    ratio = minute/60 + second/3600
    r_maritime = r_maritime_1 * ratio + r_maritime_0 * (1-ratio)
    print('Interpolated ratio of maritime aerosol in maritime/continental mixture: ' + str(r_maritime[0]))
    print('Angstrom exponent: {:.2f}'.format(r_maritime[1]))
    print('Single scattering albedo: {:.2f}'.format(r_maritime[2]))
    return r_maritime
    
    
# Find value of ancillary data at lat and lon
def _get_AER(file, lat, lon):
    
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
    
    # return ratio of maritime, rest being continental 
    def categorize_aer(i,maritime,continental):
        if maritime < continental: 
            if i < maritime: return 1
            elif i > continental: return 0
            else: return abs(i - maritime) / abs(continental - maritime)
        else:
            if i > maritime: return 1
            elif i < continental: return 0
            else: return abs(maritime - i) / abs(maritime - continental)
    
    nc = nc4.Dataset(file, 'r') 
    
    lons = linspace(nc.getncattr('geospatial_lon_min'),
                    nc.getncattr('geospatial_lon_max'),
                    nc.getncattr('geospatial_lon_resolution'))
    
    lats = linspace(nc.getncattr('geospatial_lat_min'),
                    nc.getncattr('geospatial_lat_max'),
                    nc.getncattr('geospatial_lat_resolution'))
    
    aer_ANG = nc['TOTANGSTR'][:]
    aer_SSA = nc['TOTSCATAU'][:] / nc['TOTEXTTAU'][:]
    aer_AOT = nc['TOTSCATAU'][:]
    
    interp_ANG = interpolate.RegularGridInterpolator((lats,lons), aer_ANG)
    i_ANG = interp_ANG([lat,lon])[0]
    
    interp_SSA = interpolate.RegularGridInterpolator((lats,lons), aer_SSA)
    i_SSA = interp_SSA([lat,lon])[0]
    
    interp_AOT = interpolate.RegularGridInterpolator((lats,lons), aer_AOT)
    i_AOT = interp_AOT([lat,lon])[0]
    
    # print('Angstrom exponent: ' + str(i_ANG))
    # print('SSA: ' + str(i_SSA))
    
    # according the the two 
    maritime_ANG = categorize_aer(i_ANG, maritime=0.265, continental=1.132)
    maritime_SSA = categorize_aer(i_SSA, maritime=0.98903, continental=0.89319)
    ratio_maritime = (maritime_ANG + maritime_SSA) / 2
    
    # print('ratio_maritime: ' + str(ratio_maritime))
    return np.array([ratio_maritime, i_ANG, i_SSA, i_AOT])
