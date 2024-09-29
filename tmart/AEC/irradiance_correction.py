# This file is part of T-Mart.
#
# Copyright 2024 Yulun Wu.
#
# T-Mart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# Correct for the non-linearity of the ratio of environmental irradiance to surface reflectance for homogeneous Lambertian surfaces. 

from Py6S import *

def irradiance_correction(image, wl_RC, band = None,
                          tm_vza = 0, tm_vaa = 0, tm_sza = 0, tm_saa = 0, 
                          atm_profile = None, 
                          aerosol_type = 'Maritime', aot550 = 0):
    '''
    

    Parameters
    ----------
    image : numpy array 
        surface reflectance.
    wl_RC : TYPE
        DESCRIPTION.
    band : TYPE, optional
        DESCRIPTION. The default is None.
    tm_vza : TYPE, optional
        DESCRIPTION. The default is 0.
    tm_vaa : TYPE, optional
        DESCRIPTION. The default is 0.
    tm_sza : TYPE, optional
        DESCRIPTION. The default is 0.
    tm_saa : TYPE, optional
        DESCRIPTION. The default is 0.
    atm_profile : TYPE, optional
        DESCRIPTION. The default is None.
    aerosol_type : TYPE, optional
        DESCRIPTION. The default is 'Maritime'.
    aot550 : TYPE, optional
        DESCRIPTION. The default is 0.

    Returns
    -------
    image_out : TYPE
        DESCRIPTION.

    '''
    
    import numpy as np
    import sys
    
    # An unknown bug of 6S unable to handle 551nm
    if wl_RC == 0.551: wl_RC = 0.55
    
    median_image = np.median(image)
    
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
        s.atmos_profile = AtmosProfile.UserWaterAndOzone(atm_profile['water_vapour']/10, atm_profile['ozone']/1000)

    # Aerosol 
    if aerosol_type == 'Maritime':
        aerosol_profile = s.aero_profile = AeroProfile.PredefinedType(AeroProfile.Maritime)
    elif aerosol_type == 'Continental':
        aerosol_profile = s.aero_profile = AeroProfile.PredefinedType(AeroProfile.Continental)
    elif isinstance(aerosol_type, float):
        r_mar = aerosol_type
        r_con = 1 - r_mar
        s.aero_profile = AeroProfile.User(soot = 0.01*r_con, water = 0.05*r_mar + 0.29*r_con, 
                                          oceanic = 0.95*r_mar, dust = 0.7*r_con)
    else:
        sys.exit('Unrecgnized aerosol type')
        
    s.aot550 = aot550
    
    # Elevation 
    s.altitudes.set_sensor_satellite_level()
    s.altitudes.set_target_custom_altitude(0)
    
    # Wavelength and band 
    if band is None:  
        s.wavelength = Wavelength(wl_RC) 
    else: 
        s.wavelength = band
        
    # Lists of TOA direct reflectance and corresponding surface reflectance 
    list_R_dir = []
    list_SR = [0.25,0.5,0.75,1]
    
    for SR in list_SR:
        s.ground_reflectance = GroundReflectance.HomogeneousLambertian(SR)
        s.run()
        list_R_dir.append(s.outputs.values['pixel_reflectance'])
    
    
    # Calculate the ratio 
    array_R_dir = np.array(list_R_dir)
    array_SR = np.array(list_SR)
    array_ratio = array_R_dir / array_SR
    fit_ratio_SR = np.polyfit(array_SR , array_ratio, 2)
    p = np.poly1d(fit_ratio_SR)
    median_ratio = p(median_image)
    
    # Ratio of all pixels 
    ratio = p(image)
    
    # Relative to median 
    correction = ratio / median_ratio
    
    # Do not correct pixels brighter than mean 
    correction[correction>1] = 1
    image_out = image * correction
    max_correction = correction.min()
    max_correction_percent = str(round((1 - max_correction)*100, 2))
    print('Maximum change in pixel value: ' + max_correction_percent + '%')
    
    return image_out
