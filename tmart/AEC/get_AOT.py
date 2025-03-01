# This file is part of T-Mart.
#
# Copyright 2024 Yulun Wu.
#
# T-Mart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# Calculate AOT 550 value by identifying zero-water-leaving-reflectance pixels 

def get_AOT(metadata, config, anci, mask_cloud, mask_all, n_photon):
    
    import tmart, Py6S
    import numpy as np
    import rasterio
    import time 
    from scipy import signal
    from scipy.interpolate import interp1d
    
    sensor = metadata['sensor']
    
    # Wavelength and band
    if sensor == 'S2A': 
        band = Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_8A)
        wl = 864.7
        band_file = metadata['B8A']
    elif sensor == 'S2B': 
        band = Py6S.Wavelength(Py6S.PredefinedWavelengths.S2B_MSI_8A)
        wl = 864.0
        band_file = metadata['B8A']
    elif sensor == 'S2C':
        band = Py6S.Wavelength(0.83, 0.9, [0.0, 8.291e-05, 0.00016582, 0.00018964, 0.00021346, 0.00032219, 0.00043092, 0.002265495, 0.00410007, 0.15940945, 0.31471883, 0.6560948, 0.99747077, 0.980849005, 0.96422724, 0.950272335, 0.93631743, 0.7946804599999999, 0.65304349, 0.33604138499999997, 0.01903928, 0.010021164999999999, 0.00100305, 0.0006311349999999999, 0.00025922, 0.000206555, 0.00015389, 7.6945e-05, 0.0])
        wl = 865.6
        band_file = metadata['B8A']
    elif sensor == 'L8':
        band = Py6S.Wavelength(Py6S.PredefinedWavelengths.LANDSAT_OLI_B5)
        wl = 864.7
        band_file = metadata['B5']
    
    # Load files
    band_ds = rasterio.open(band_file)
    
    # Resolution of this particular band 
    res_band = int(abs(band_ds.transform[0]))
    res_AEC = int( metadata['resolution'] * metadata['reshape_factor'])
    
    # Full resolution to AEC resolution 
    height_reshaped = int( metadata['AEC_height'] / metadata['reshape_factor']) 
    width_reshaped = int( metadata['AEC_width'] / metadata['reshape_factor']) 

    # Number of columns and rows to pad
    pad_columns = metadata['AEC_width'] - metadata['width'] 
    pad_rows = metadata['AEC_height'] - metadata['height'] 
    pad_rows_tmp = int(pad_rows/(res_band/metadata['resolution']))
    pad_columns_tmp = int(pad_columns/(res_band/metadata['resolution']))
    
    # Scaling  
    if sensor == 'S2A' or sensor == 'S2B':
        scale_mult = metadata['B8A_mult']
        scale_add = metadata['B8A_add']
    elif sensor == 'L8': 
        scale_mult = metadata['B5_mult']
        scale_add = metadata['B5_add']
    
    # Read image
    image = band_ds.read(1)
    image = np.pad(image, ((0, pad_rows_tmp), (0, pad_columns_tmp)), mode='constant', constant_values=0)
    is_nan = image==0
    image = image * scale_mult + scale_add
    
    # Turn negative to 0 
    image[image<0] = 0
    
    # Mask no_data
    image[is_nan] = np.nan
    
    # Fill nan with closest values
    image = tmart.AEC.fillnan(image)
    
    # Reshape_factor specific to this resolution 
    reshape_factor_tmp = int(metadata['reshape_factor'] * metadata['resolution']/res_band)
    
    # Reshape
    image_AEC = image.reshape([height_reshaped, reshape_factor_tmp, 
                               width_reshaped, reshape_factor_tmp]).mean(3).mean(1)
    
    # Loop through AOTs
    AOTs = [0.0, 0.1, 0.3]

    # Record lowest R at each of the AOTs
    R_lowest = []
    
    for AOT in AOTs:
    
        AEC_parameters = tmart.AEC.get_parameters(n_photon = n_photon, SR = np.nanmedian(image), wl = wl, band = band, 
                                                  target_pt_direction=metadata['tm_pt_dir'], sun_dir=metadata['tm_sun_dir'], 
                                                  atm_profile = anci, 
                                                  aerosol_type = anci['r_maritime'], aot550 = AOT, 
                                                  cell_size = res_AEC,
                                                  window_size = metadata['window_size'], isWater = 1)
    
        conv_window_1   = AEC_parameters['conv_window_1']
        F_correction    = AEC_parameters['F_correction']
        F_captured      = AEC_parameters['F_captured']
        R_atm           = AEC_parameters['R_atm']
        R_glint         = AEC_parameters['R_glint']
        
        # Cap glint at max_glint
        if R_glint > float(config['max_glint']): R_glint = float(config['max_glint'])
        
        # image_R_surf: reshaped 
        image_R_surf = image_AEC - R_atm
        print('\nNumber of negative R_surf pixels: ' + str(np.sum(image_R_surf < 0)) + '/' + str(height_reshaped * width_reshaped))
        
        # Turn negative to 0
        image_R_surf[image_R_surf<0] = 0
    
        # Read AEC parameters 
        hp_cloud_contribution = float(config['cloud_contribution'])
        hp_AE_land = config['AE_land'] == 'True'
        
        # Reduce contribution of clouds according to setting 
        image_R_surf[mask_cloud[str(res_AEC) + 'm']] = image_R_surf[mask_cloud[str(res_AEC) + 'm']] * hp_cloud_contribution
        
        # Convolution 
        start_time = time.time()
        print("\nConvolution started ")
        filter_kernel = np.flip(conv_window_1) # it's flipped in convolve by default
        R_conv = signal.convolve2d(image_R_surf, filter_kernel,
                                         mode='same', boundary='fill', fillvalue=image_R_surf.mean())
        print("Convolution completed: %s seconds " % (time.time() - start_time))
        
        # Smoothing the edges 
        if reshape_factor_tmp>1 and not hp_AE_land: 
            image_water = np.where(mask_all[str(res_band) + 'm'], np.nan, image)
            
            # Suppress warning of mean of empty slice 
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                image_water_AEC = image_water.reshape([height_reshaped, reshape_factor_tmp, width_reshaped, reshape_factor_tmp])
                image_water_AEC = np.nanmean(image_water_AEC,3)
                image_water_AEC = np.nanmean(image_water_AEC,1)
    
            # with non-water as nan 
            image_water_R_surf = image_water_AEC - R_atm
            
            # Edit image_R_surf: where there's water image_R_surf, keep it, where there's nan, use original image_R_surf
            # image_R_surf can only be edited after convolution 
            image_R_surf = np.where(np.isnan(image_water_R_surf), image_R_surf, image_water_R_surf)   
    
        R_correction = (R_conv - image_R_surf) * F_correction 
        print('\nNumber of pixels where R_correction > R_surf : ' + str(np.sum(R_correction>image_R_surf)) + '/' + str(height_reshaped * width_reshaped))
        
        # Back to the original size 
        R_correction_original_shape = np.repeat(np.repeat(R_correction, reshape_factor_tmp, axis=0), 
                                                reshape_factor_tmp, axis=1)
        
        # Skipping smoothing the gridline artifacts 
        
        # User specified minimum R of water in NIR 
        min_R_NIR = float(config['min_R_NIR'])
        
        # remove AE from TOA reflectance 
        temp_SR = image - R_correction_original_shape - R_atm
        
        # Correcting for non-linearity of the ratio of environmental irradiance to surface reflectance for homogeneous Lambertian surfaces
        # This is implemented for water only
        # Input has to include land, because the function uses the env irradiance of the average reflectance across the scene
        print("\nIrradiance correction: ")
        temp_SR_water = tmart.AEC.irradiance_correction(image = temp_SR, wl_RC = wl/1000, band = band,
                                                        tm_vza = metadata['vza'], tm_vaa = metadata['vaa'], 
                                                        tm_sza = metadata['sza'], tm_saa = metadata['saa'],
                                                        atm_profile = anci, 
                                                        aerosol_type = anci['r_maritime'], aot550 = AOT)
        
        # Residual water-leaving reflectance, we want the darkest pixels to be near 0
        temp_SR_water = temp_SR_water - R_glint - min_R_NIR
        
        # Add nan back 
        temp_SR_water[is_nan] = np.nan
        
        # Flatten the 2D array to a 1D array for plotting
        flattened_data = temp_SR_water.flatten()
        
        sorted_data = np.sort(flattened_data)
        smallest_values = sorted_data[:10_000]
        
        
        
        
        ### Consider an adaptive number instead of 10_000 
        # Either by area or correlation coefficient
        
        
        
        # Create an array of indices for the smallest values
        indices = np.arange(10_000)
        
        
        '''
        ### Scatter plot 
        
        from matplotlib import pyplot as plt
        
        # Create a scatter plot
        plt.scatter(indices, smallest_values, marker='o', s=10)  # Adjust marker and size as needed
        plt.xlabel('Index')
        plt.ylabel('Value')
        plt.title('Scatter Plot of Smallest Values, AOT550: ' + str(AOT))
        plt.grid(True)
        
        # Show the plot
        plt.show()
        '''
            
    
        # Regression 
        from scipy import stats
        slope, intercept, r_value, p_value, std_err = stats.linregress(indices, smallest_values)
        print('\nLowest reflectance: ' + str(intercept))
        R_lowest.append(intercept)
    
    print('\nAOTs: ' + str(AOTs))
    print('R_lowest: ' + str(R_lowest))
    
    # Interpolate AOT
    if min(R_lowest) > 0: interpolated_AOT = max(AOTs)
    elif max(R_lowest) < 0: interpolated_AOT = min(AOTs)
    else:
        # prioritize the first two values, in case the line crosses 0 twice 
        if R_lowest[1] < 0:    
            AOTs = AOTs[ : -1]
            R_lowest = R_lowest[ : -1]
        
        # Create an interpolation function
        f = interp1d(R_lowest, AOTs, kind='linear')
        
        # Find x when y is 0
        interpolated_AOT = f(0)

    # Limit the range 
    if interpolated_AOT > max(AOTs): interpolated_AOT = max(AOTs)
    if interpolated_AOT < min(AOTs): interpolated_AOT = min(AOTs)
    print('interpolated_AOT: ' + str(interpolated_AOT))
    return interpolated_AOT 


