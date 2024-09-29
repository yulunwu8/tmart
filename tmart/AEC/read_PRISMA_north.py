# This file is part of T-Mart.
#
# Copyright 2024 Yulun Wu.
#
# T-Mart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# read PRSIMA north - return the top direction relative to north 

def read_PRISMA_north(file): 
    # assume PRISMA images are 1000x1000

    import pyproj
    import netCDF4 as nc4
    
    with nc4.Dataset(file, 'r')  as dset:
    
        a_lat = dset['lat'][:]
        a_lon = dset['lon'][:]
        
        # image = dset['rhot_505'][:]
    
    # lat lon in four corners 
    top_left = [a_lat[0,0],a_lon[0,0]]
    top_right = [a_lat[0,999],a_lon[0,999]]
    bot_left = [a_lat[999,0],a_lon[999,0]]
    bot_right = [a_lat[999,999],a_lon[999,999]]
    
    geodesic = pyproj.Geod(ellps='WGS84')
    
    up_left,__,__ = geodesic.inv(bot_left[1], bot_left[0], top_left[1], top_left[0])
    # print(up_left)
    
    up_right,__,__ = geodesic.inv(bot_right[1], bot_right[0], top_right[1], top_right[0])
    # print(up_right)
    
    return (up_left + up_right) / 2

