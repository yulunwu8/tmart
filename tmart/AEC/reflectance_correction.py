# This file is part of T-Mart.
#
# Copyright 2023 Yulun Wu.
#
# T-Mart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


### Correct for the non-linearity of the ratio of environmental irradiance to surface reflectance for homogeneous Lambertian surfaces. 

from Py6S import *

def reflectance_correction(image, wl_RC, 
                           tm_vza = 0, tm_vaa = 0, tm_sza = 0, tm_saa = 0, 
                           atm_profile = None, 
                           aerosol_type = 'Maritime', aot550 = 0):
    
    import numpy as np
    import sys
    
    
    # An unknown bug of 6S unable to handle 551nm
    if wl_RC == 0.551: wl_RC = 0.55
    
    
    mean_image = np.mean(image)
    
    # SixS
    s = SixS()
    
    s.geometry = Geometry.User()
    s.geometry.solar_z = tm_sza
    s.geometry.solar_a = tm_saa
    s.geometry.view_z = tm_vza
    s.geometry.view_a = tm_vaa
    
    # Atmosphere
    if atm_profile is None:
        s.atmos_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeSummer)
    else: 
        s.atmos_profile = AtmosProfile.UserWaterAndOzone(atm_profile[0], atm_profile[1])
     
    # Aerosol 
    if aerosol_type == 'Maritime':
        aerosol_profile = s.aero_profile = AeroProfile.PredefinedType(AeroProfile.Maritime)
    elif aerosol_type == 'Continental':
        aerosol_profile = s.aero_profile = AeroProfile.PredefinedType(AeroProfile.Continental)
    else:
        sys.exit('Unrecgnized aerosol type')
        
    s.aot550 = aot550
    
    s.altitudes.set_sensor_satellite_level()
    s.altitudes.set_target_custom_altitude(0)
    
    
    s.wavelength = Wavelength(wl_RC) 
    
    # Lists of TOA direct reflectance and corresponding surface reflectance 
    list_R_dir = []
    list_SR = [0.25,0.5,0.75,1]
    
    for SR in list_SR:
        s.ground_reflectance = GroundReflectance.HomogeneousLambertian(SR)
        s.run()
        list_R_dir.append(s.outputs.values['pixel_reflectance'])
    
    list_R_dir
    
    array_R_dir = np.array(list_R_dir)
    array_SR = np.array(list_SR)
    array_ratio = array_R_dir / array_SR
    
    fit_ratio_SR = np.polyfit(array_SR , array_ratio, 2)
    p = np.poly1d(fit_ratio_SR)
    mean_ratio = p(mean_image)
    ratio = p(image)
    
    correction = ratio / mean_ratio
    
    # do not correct pixels brighter than mean 
    correction[correction>1] = 1
    
    image_out = image * correction
    
    max_correction = correction.min()
    
    # max_correction = np.nanmin(correction)
    
    max_correction_percent = str(round((1 - max_correction)*100, 2))
    
    print('Maximum change: ' + max_correction_percent + '%')
    
    
    return image_out




# if __name__ == "__main__":
#     import numpy as np
#     import netCDF4 as nc4
#     file = 'test.nc'
#     with nc4.Dataset(file, 'r')  as dset:
#         a_test = dset['rhot_551'][:]
#     test = reflectance_correction(a_test, 0.551)
    
    

    





















