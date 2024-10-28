# This file is part of T-Mart.
#
# Copyright 2024 Yulun Wu.
#
# T-Mart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# Derive AE correction parameters 

def get_parameters(n_photon = 10_000, SR = 0.5, 
                   wl = 833, band = None, 
                   target_pt_direction=[180,0], sun_dir=[0,0], 
                   atm_profile = None, 
                   aerosol_type = 'Maritime', aot550 = 0.2, 
                   cell_size = 100,window_size = None,
                   window_size_x = None, window_size_y = None, isWater = 0,
                   njobs=100):
    
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
        atm_profile = AtmosProfile.UserWaterAndOzone(atm_profile['water_vapour']/10, atm_profile['ozone']/1000)
        
    # DEM and reflectance 
    image_DEM = np.full((window_size_x, window_size_y), 0) # in meters
    image_reflectance = np.full((window_size_x, window_size_y), SR) # unitless     
    image_isWater = np.full((window_size_x, window_size_y), isWater) # 1 is water, 0 is land
    
    # Synthesize a surface object
    my_surface = tmart.Surface(DEM = image_DEM,
                               reflectance = image_reflectance,
                               isWater = image_isWater,
                               cell_size = cell_size)  
    my_surface.set_background(bg_isWater=isWater)                               
    
    ### Atmosphere ###
    my_atm = tmart.Atmosphere(atm_profile, aot550 = aot550, aerosol_type = aerosol_type)
    
    ### Running T-Mart ###
    my_tmart = tmart.Tmart(Surface = my_surface, Atmosphere= my_atm, shadow=False)
    my_tmart.set_wind(wind_speed=1, wind_azi_avg = True)
    
    my_tmart.set_geometry(target_pt_direction=target_pt_direction,
                          pixel=[int(window_size_y/2),int(window_size_x/2)], 
                          sun_dir=sun_dir)    
    
    results = my_tmart.run(wl=wl, band=band, n_photon=n_photon, njobs=njobs)
    # results = my_tmart.run_plot(wl=wl, plot_on=True, plot_range=[0,cell_size*window_size_x,0,cell_size*window_size_x,0,100_000])
    
    # Calculate reflectances using recorded photon information 
    R = tmart.calc_ref(results,detail=True)
    for k, v in R.items():
        print(k, '     ' , v)
        
    ### Computing parameters  
    
    # column names 
    columns = ['pt_id', 'movement', 'L_cox-munk', 'L_whitecap', 'L_water', 'L_land', 
               'L_rayleigh', 'L_mie', 'x', 'y', 'z', 'shadowed', 'if_env']
    
    # Action: make this a numpy array to speed up computation 
    df = pd.DataFrame(results, columns = columns)
    df_env = df[df.if_env==1].copy()
    
    df_env_sum = df_env.iloc[:,2:6].sum(axis=1)
    df_env['L_surface'] = df_env_sum
    
    ### Bin points to convolution matrix  
    
    # Old method 
    
    # image_env = np.full((window_size_x, window_size_y), 0.0)
    # for xi in range(window_size_x):
    #     for yi in range(window_size_y):
    #         print('xi: ' + str(xi))
    #         print('yi: ' + str(yi))
    #         xmin = xi * cell_size 
    #         xmax = (xi + 1) * cell_size 
    #         ymin = yi * cell_size
    #         ymax = (yi + 1) * cell_size
    #         r_sum = df_env[(df_env.x>xmin) & (df_env.y>ymin) & 
    #                         (df_env.x<xmax) & (df_env.y<ymax)].L_surface.sum() # .L_land.sum()
    #         image_env[yi,xi] = r_sum 
    
    # Classify df_env rows into cells
    x_bins = np.linspace(0, cell_size * window_size_x, window_size_x + 1)
    y_bins = np.linspace(0, cell_size * window_size_y, window_size_y + 1)

    # Use np.histogram2d to compute the sum of values in each cell
    image_env, _, _ = np.histogram2d(df_env['y'], df_env['x'], bins=[y_bins, x_bins], weights=df_env['L_surface'])
    
    conv_window = image_env.copy()
    
    F_captured = conv_window.sum()/n_photon / R['R_env'] # ratio of R_env included in image_env
    
    print('\nR_env captured in conv_window: ' + str(F_captured))

    ### normalize the sum of the remaining pixels to 1 
    conv_window_1 = conv_window.copy()
    # multiply centre pixel by F_captured
    conv_window_1[int(conv_window_1.shape[0]/2),
                  int(conv_window_1.shape[1]/2)] = F_captured * conv_window_1[int(conv_window_1.shape[0]/2),
                                                                              int(conv_window_1.shape[1]/2)]
    conv_window_1 = conv_window_1 / conv_window_1.sum()
    
    ### correction factor 
    F_correction = (R['R_env'] / R['R_dir'])  * (1 - conv_window_1[int(conv_window_1.shape[0]/2),int(conv_window_1.shape[1]/2)])
    print('alpha: ' + str(F_correction)) # alpha is the term used in Wu et al. 2024
    
    return {'conv_window_1': conv_window_1,
            'F_correction': F_correction,
            'F_captured': F_captured,
            'R_atm': R['R_atm'],
            'R_glint': R['_R_dir_coxmunk'] + R['_R_env_coxmunk']  }
