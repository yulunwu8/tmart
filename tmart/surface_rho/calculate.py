# This file is part of T-Mart.
#
# Copyright 2023 Yulun Wu.
#
# T-Mart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

def calculate(wl, viewing_zenith, solar_zenith, relative_azimuth, 
              aot550=0.05, wind_speed=3,n_photon = 100_000,
              atm_profile=None,aerosol_type=None,spectral_surface=None,as_pandas_df=True):
    
    '''Calculate sea-surface reflectance factor. See Wu et al. 2023, Lee et al. 2010, or Mobley 1999 for full description. 

    Arguments:
    
    * ``wl`` -- central wavelength in nm, or a list for iteration: [starting_wavelength, ending wavelength, interval].
    * ``viewing_zenith`` -- viewing zenith angle in degrees.
    * ``solar_zenith`` -- solar zenith angle in degrees.
    * ``relative_azimuth`` -- relative azimuth angle in degrees.
    * ``aot550`` -- aerosol optical thickness at 550nm.
    * ``wind_speed`` -- wind speed in m/s.
    * ``n_photon`` -- number of photons in each T-Mart run. 100_000 is recommended, 10_000 can be used to view quick results. 
    * ``atm_profile`` -- AtmosProfile object from Py6S. Default 'MidlatitudeSummer'
    * ``aerosol_type`` -- 'BiomassBurning', 'Continental', 'Desert', 'Maritime', 'Stratospheric' or 'Urban', as provided by 6S. Default 'Maritime'
    * ``spectral_surface`` -- a SpectralSurface object in T-Mart. Default 'water_chl1'
    * ``as_pandas_df`` -- If true, return a pandas dataframe, else return a dictionary
    
    Example usage:: 
    
      tmart.surface_r.calculate(wl=800, viewing_zenith=40, solar_zenith=30, relative_azimuth=135)

    '''
    
    import tmart
    import numpy as np
    import pandas as pd
    from Py6S.Params.atmosprofile import AtmosProfile
    
    ### Input 
    
    if atm_profile==None: 
        atm_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeSummer) 
    
    if aerosol_type==None: 
        aerosol_type = 'Maritime' 
    
    if spectral_surface == None:
        spectral_surface = tmart.SpectralSurface('water_chl1')
    
    rho_output = []
    
    if isinstance(wl, list):
        wls = range(wl[0],wl[1]+wl[2],wl[2])
    else:
        wls = [wl]
    
    for wl in wls:
    
        ### L_sky ###
        
        # Three same-size numpy arrays are needed
        image_DEM = np.array([[0,0],[0,0]]) # in meters
       
        image_reflectance = np.array([[0,0],[0,0]]) # unitless     
        
        image_isWater = np.array([[1,1],[1,1]])# 1 is water, 0 is land
        
        # pixel width and length, in meters 
        cell_size = 0.1
        
        # Synthesize a surface object
        my_surface = tmart.Surface(image_DEM,image_reflectance,image_isWater,cell_size)  
        # Set background information, 1 or 2 background surfaces can be set;
        # If 2: the first background is the one closer to [0,0]
        my_surface.set_background(bg_ref        = spectral_surface.wl(wl), 
                                  bg_isWater    = 1,
                                  bg_elevation  = 0, # elevation of both background
                                  bg_coords     = [[0,0],[10,10]]) # a line dividing two background
        
        # Synthesize an atmosphere object     
        my_atm = tmart.Atmosphere(atm_profile, aot550, aerosol_type)
        
        ### Running T-Mart ###
        my_tmart = tmart.Tmart(Surface = my_surface, Atmosphere= my_atm)
        my_tmart.set_wind(wind_speed=wind_speed,wind_azi_avg=True)
        my_tmart.set_water(water_salinity=35, water_temperature=20)
        my_tmart.set_geometry(sensor_coords=[0.075,0.05,0.001], 
                              target_pt_direction=[viewing_zenith,relative_azimuth],
                              sun_dir=[solar_zenith,0])
        
        results = my_tmart.run(wl=wl, n_photon=n_photon)
        R1 = tmart.calc_ref(results, n_photon=n_photon)
        R1['Wavelength'] = wl
        R1['Type'] = 'L_sky'
        for k, v in R1.items():
    	    print(k, '     ' , v)
        
        ### L_sr ###
        
        my_tmart.set_geometry(sensor_coords=[0.075,0.05,0.001], 
                              target_pt_direction=[180-viewing_zenith,relative_azimuth],
                              sun_dir=[solar_zenith,0])
        
        
        results = my_tmart.run(wl=wl, n_photon=n_photon)
        R2 = tmart.calc_ref(results, n_photon=n_photon)
        R2['Wavelength'] = wl
        R2['Type'] = 'L_sr'
        for k, v in R2.items():
    	    print(k, '     ' , v)
        
        
        ### calculate rho ###
        
        rho = R2['R_total'] / R1['R_total']
        print('Sea-surface reflectance factor: ' + str(rho))
        rho_output.append({'wavelength': wl, 'rho': rho})
    
    if as_pandas_df: rho_output = pd.DataFrame(rho_output)
    
    return rho_output 




