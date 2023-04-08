# This file is part of TMart.
#
# Copyright 2023 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.




from Py6S import *


def reflectance_correction(image, wl_RC, 
                           tm_vza, tm_vaa, tm_sza, tm_saa, 
                           atm_profile = None, 
                           aerosol_type = 'Maritime', aot550 = 0.2):
    
    import numpy as np
    import sys
    
    

    mean_image = np.mean(image)
    
        
    # SixS

    s = SixS()
    
    s.geometry = Geometry.User()
    s.geometry.solar_z = tm_sza
    s.geometry.solar_a = tm_saa
    
    s.geometry.view_z = tm_vza
    s.geometry.view_a = tm_vaa
    
    if atm_profile is None:
        atm_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeSummer)
    else: 
        atm_profile = AtmosProfile.UserWaterAndOzone(atm_profile[0], atm_profile[1])
        
    s.atmos_profile = atm_profile
     
    
    if aerosol_type == 'Maritime':
        aerosol_profile = s.aero_profile = AeroProfile.PredefinedType(AeroProfile.Maritime)
    elif aeaerosol_type == 'Continental':
        aerosol_profile = s.aero_profile = AeroProfile.PredefinedType(AeroProfile.Continental)
    else:
        sys.exit('Unrecgnized aerosol type')
        
    s.aot550 = aot550
    
    s.altitudes.set_sensor_satellite_level()
    s.altitudes.set_target_custom_altitude(0)
    
    
    s.wavelength = Wavelength(wl_RC) 
    
    list_R_dir = []
    
    list_SR = [0.25,0.5,0.75,1]
    
    for SR in list_SR:
        s.ground_reflectance = GroundReflectance.HomogeneousLambertian(SR)
        s.run()
        
        # print('\n===background_reflectance:')
        # print (s.outputs.values['background_reflectance'])
        
        # print('\n===pixel_reflectance:')
        # print (s.outputs.values['pixel_reflectance'])
        
        # print('\n===atmospheric_intrinsic_reflectance:')
        # print (s.outputs.values['atmospheric_intrinsic_reflectance'])
        
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
    
    
    # print('\n===image_out:')
    # print(image_out)
    
    # print('\n===correction:')
    # print(correction)
    
    max_correction = correction.min()
    print('Reflectance correction, maximum: ' + str(max_correction))
    
    
    return image_out


if __name__== "__main__" :
    
    import numpy as np

    image = np.arange(16).reshape(4,4) / 32.
    
    print('\n===image:')
    
    print(image)
    
    
    
    # mean_image = np.mean(image)
    # print('mean: ' + str(mean_image))

    
    tm_vza = 0
    tm_vaa = 0
    tm_sza = 0
    tm_saa = 0
    
    test = reflectance_correction(image = image, wl_RC = 0.4, tm_vza = tm_vza, tm_vaa = tm_vaa, tm_sza = tm_sza, tm_saa = tm_saa)





































