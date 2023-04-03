# This file is part of TMart.
#
# Copyright 2023 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


### Derive AE correction parameters 

def get_parameters(n_photon = 10_000, SR_avg = 1, wl = 833, band = None, 
        target_pt_direction=[180,0], sun_dir=[0,0], 
        atm_profile = None, 
        aerosol_type = 'Maritime', aot550 = 0.2, 
        cell_size = 100,window_size = None,
        window_size_x = None, window_size_y = None):
    
    
    import tmart
    import numpy as np
    import pandas as pd 
    from Py6S.Params.atmosprofile import AtmosProfile

    if window_size is not None: 
        window_size_x = window_size
        window_size_y = window_size
           
    ### Action: add test 
    # window_size_y and window_size_y have to be odd integers  
    
    if atm_profile is None:
        atm_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeSummer)
    else: 
        atm_profile = AtmosProfile.UserWaterAndOzone(atm_profile[0], atm_profile[1])
        
    # DEM and reflectance 
    image_DEM = np.full((window_size_x, window_size_y), 0) # in meters
    image_reflectance = np.full((window_size_x, window_size_y), SR_avg) # unitless     
    image_isWater = np.full((window_size_x, window_size_y), 0) # 1 is water, 0 is land
    
    # Synthesize a surface object
    my_surface = tmart.Surface(DEM = image_DEM,
                               reflectance = image_reflectance,
                               isWater = image_isWater,
                               cell_size = cell_size)  
                                   
    ### Atmosphere ###
    my_atm = tmart.Atmosphere(atm_profile, aot550 = aot550, aerosol_type = aerosol_type)
    
    ### Running T-Mart ###
    my_tmart = tmart.Tmart(Surface = my_surface, Atmosphere= my_atm, shadow=False)
    
    my_tmart.set_geometry(target_pt_direction=target_pt_direction,
                          pixel=[int(window_size_y/2),int(window_size_x/2)], 
                          sun_dir=sun_dir)    
    
    results = my_tmart.run(wl=wl, band=band, n_photon=n_photon)
    # results = my_tmart.run_plot(wl=wl, plot_on=True, plot_range=[0,cell_size*window_size_x,0,cell_size*window_size_x,0,100_000])
    
    # Calculate reflectances using recorded photon information 
    R = tmart.calc_ref(results)
    for k, v in R.items():
        print(k, '     ' , v)
        
    ### Computing parameters  
    
    # column names 
    columns = ['pt_id', 'movement', 'L_cox-munk', 'L_whitecap', 'L_water', 'L_land', 
               'L_rayleigh', 'L_mie', 'x', 'y', 'z', 'shadowed', 'if_env']
    
    # Action: make this a numpy array to speed up computation 
    df = pd.DataFrame(results, columns = columns)
    df_env = df[df.if_env==1]
    
    ### Bin points to convolution matrix  
    image_env = np.full((window_size_x, window_size_y), 0.0)
    for xi in range(window_size_x):
        for yi in range(window_size_y):
            xmin = xi * cell_size 
            xmax = (xi + 1) * cell_size 
            ymin = yi * cell_size
            ymax = (yi + 1) * cell_size
            r_sum = df_env[(df_env.x>xmin) & (df_env.y>ymin) & 
                           (df_env.x<xmax) & (df_env.y<ymax)].L_land.sum()
            image_env[yi,xi] = r_sum 
    
    conv_window = image_env.copy()
    print('\nR_env captured in conv_window: ')
    F_captured = conv_window.sum()/n_photon / R['R_env'] # ratio of R_env included in image_env
    print(F_captured)
    
    ### normalize the sum of the remaining pixels to 1 
    conv_window_1 = conv_window.copy()
    # multiply centre pixel by F_captured
    conv_window_1[int(conv_window_1.shape[0]/2),
                  int(conv_window_1.shape[1]/2)] = F_captured * conv_window_1[int(conv_window_1.shape[0]/2),
                                                                              int(conv_window_1.shape[1]/2)]
    conv_window_1 = conv_window_1 / conv_window_1.sum()
    
    ### correction factor 
    F_correction = R['R_env'] / R['R_dir']
    print('\nF_correction: ')
    print(F_correction)
    
    return {'conv_window_1': conv_window_1,
            'F_correction': F_correction,
            'F_captured': F_captured,
            'R_atm': R['R_atm']}




'''

if __name__ == "__main__":
   
    import Py6S
    import tmart

    wl = 833
    band = Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_08)
    aerosol_type = 'Continental'
    aot550 = 0.28
    
    test = tmart.AEC.get_parameters(wl = wl, band = band,
                                    aerosol_type = aerosol_type, aot550 = aot550,
                                    cell_size = 200,
                                    window_size=5)
    
    print(test)
    test2 = test['conv_window_1']
    

### Test AEC

# image_test = image_reflectance.copy()

image_test = np.full((window_size_x, window_size_y), 0.39378) # unitless    


# image_test[:,6:] = 0.5

image_test[5,5] = 0

import scipy.signal

filter_kernel = np.flip(conv_window_1) # it's flipped in convolve by default
R_conv = scipy.signal.convolve2d(image_test, filter_kernel,
                              mode='same', boundary='fill', fillvalue=image_test.mean())


# print(R_conv)

# centre 
centre_R_conv = R_conv[int(R_conv.shape[0]/2),int(R_conv.shape[1]/2)]
print('\ncentre_R_conv:')
print(centre_R_conv)


'''
























