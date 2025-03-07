# This file is part of T-Mart.
#
# Copyright 2024 Yulun Wu.
#
# T-Mart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# AEC for a single wavelength or band 

def AEC(AEC_band_name, AEC_band_6S, wl, AOT, metadata, config, anci, mask_cloud, mask_all, n_photon, njobs):
    
    import tmart
    import rasterio, sys, time, gc
    from scipy import signal
    import numpy as np
    import math
    
    sensor = metadata['sensor']
    
    # Load files 
    band_file = metadata[AEC_band_name]
    band_ds = rasterio.open(band_file, mode='r+', QUALITY='100', REVERSIBLE='yes')
    
    # Resolution of this particular band 
    res_band = int(abs(band_ds.transform[0]))
    res_AEC = int( metadata['resolution'] * metadata['reshape_factor'])
    
    # Full resolution to AEC resolution (padded resolution)
    height_reshaped = int( metadata['AEC_height'] / metadata['reshape_factor']) 
    width_reshaped = int( metadata['AEC_width'] / metadata['reshape_factor']) 

    # Number of columns and rows to pad
    pad_columns = metadata['AEC_width'] - metadata['width'] 
    pad_rows = metadata['AEC_height'] - metadata['height'] 
    pad_rows_tmp = int(pad_rows/(res_band/metadata['resolution']))
    pad_columns_tmp = int(pad_columns/(res_band/metadata['resolution']))
    
    # Scaling
    scale_mult = metadata[str(AEC_band_name) + '_mult']
    scale_add = metadata[str(AEC_band_name) + '_add']
    
    # Read image
    image = band_ds.read(1)
    image = np.pad(image, ((0, pad_rows_tmp), (0, pad_columns_tmp)), mode='constant', constant_values=0)
    is_nan = image==0
    image = image * scale_mult + scale_add
    
    # L8 solar zenith correction: https://www.usgs.gov/landsat-missions/using-usgs-landsat-level-1-data-product
    if sensor =='L8' or sensor == 'L9': image = image / math.cos(metadata['sza']/180*math.pi)
        
    # Turn negative TOA to 0 
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

    # Calculate AEC parameters 
    AEC_parameters = tmart.AEC.get_parameters(n_photon = n_photon, SR = np.nanmedian(image), wl = wl, band = AEC_band_6S, 
                                              target_pt_direction=metadata['tm_pt_dir'], sun_dir=metadata['tm_sun_dir'], 
                                              atm_profile = anci, 
                                              aerosol_type = anci['r_maritime'], aot550 = AOT, 
                                              cell_size = res_AEC,
                                              window_size = metadata['window_size'], isWater = 0, njobs=njobs)
    
    conv_window_1   = AEC_parameters['conv_window_1']
    F_correction    = AEC_parameters['F_correction']
    F_captured      = AEC_parameters['F_captured']
    R_atm           = AEC_parameters['R_atm']

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
    R_conv = signal.convolve2d(image_R_surf, filter_kernel, mode='same', boundary='fill', fillvalue=image_R_surf.mean())
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
      
    # Correction to make 
    R_correction = (R_conv - image_R_surf) * F_correction 
    print('\nNumber of pixels where R_correction > R_surf: ' + str(np.sum(R_correction>image_R_surf)) + '/' + str(height_reshaped * width_reshaped))
    
    # Back to the original size 
    R_correction_original_shape = np.repeat(np.repeat(R_correction, reshape_factor_tmp, axis=0), reshape_factor_tmp, axis=1)
    
    # Smoothing the gridline artifacts 
    if reshape_factor_tmp>1: 
    
        # Keep only water, short name for R_correction_original_shape_water
        RC_water = np.where(mask_all[str(res_band) + 'm'], np.nan, R_correction_original_shape)
        RC_water_fill_nan = tmart.AEC.fillnan(RC_water)
        
        filter_size = 3
        filter_kernel = np.full((filter_size,filter_size), 1/(filter_size**2))
        
        RC_water_smooth = signal.convolve2d(RC_water_fill_nan, filter_kernel,
                                                  mode='same', boundary='fill', fillvalue=np.nanmean(RC_water))
        
        # place RC_water_smooth in R_correction_original_shape where there's water
        R_correction_original_shape = np.where(mask_all[str(res_band) + 'm'], R_correction_original_shape, RC_water_smooth)
    
    # remove AE from TOA reflectance 
    temp_SR = image - R_correction_original_shape - R_atm
    
    # Correcting for non-linearity of the ratio of environmental irradiance to surface reflectance for homogeneous Lambertian surfaces
    # This is implemented for water only
    # Input has to include land, because the function uses the env irradiance of the average reflectance across the scene
    print("\nIrradiance correction: ")
    temp_SR_water = tmart.AEC.irradiance_correction(image = temp_SR, wl_RC = wl/1000, band = AEC_band_6S,
                                                    tm_vza = metadata['vza'], tm_vaa = metadata['vaa'], 
                                                    tm_sza = metadata['sza'], tm_saa = metadata['saa'],
                                                    atm_profile = anci, 
                                                    aerosol_type = anci['r_maritime'], aot550 = AOT)
    
    # If land correction is needed 
    if hp_AE_land:
        temp_out = np.where(mask_all[str(res_band) + 'm'], temp_SR, temp_SR_water)
        
        # Back to original dimension 
        if pad_rows_tmp>0: temp_out = temp_out[:-pad_rows_tmp]
        if pad_columns_tmp>0: temp_out = temp_out[:, :-pad_columns_tmp]
        
        # Negative to 0, this ensures the lowest TOA reflectance is no lower than R_atm
        # temp_out[temp_out<0] = 0
        
        # Scaling
        temp_out = temp_out + R_atm # TOA reflectance
        
        if sensor =='L8' or sensor == 'L9':
            temp_out = np.maximum(((temp_out * math.cos(metadata['sza']/180*math.pi) - scale_add) / scale_mult),1).astype(int)
        else:
            temp_out = np.maximum(((temp_out - scale_add) / scale_mult),1).astype(int)
    
    # If only water correction 
    else:
        # Negative to 0
        # temp_SR_water[temp_SR_water<0] = 0
        
        temp_out = temp_SR_water + R_atm # TOA reflectance
        temp_mask = mask_all[str(res_band) + 'm']
        
        # back to original dimension 
        if pad_rows_tmp>0: 
            temp_out = temp_out[:-pad_rows_tmp]
            temp_mask = temp_mask[:-pad_rows_tmp]
        if pad_columns_tmp>0: 
            temp_out = temp_out[:, :-pad_columns_tmp]
            temp_mask = temp_mask[:, :-pad_columns_tmp]
        
        # Scaling
        if sensor =='L8' or sensor == 'L9':
            temp_out = np.where(temp_mask,    band_ds.read(1),  np.maximum(((temp_out * math.cos(metadata['sza']/180*math.pi) - scale_add) / scale_mult),1).astype(int) )
        else:
            temp_out = np.where(temp_mask,    band_ds.read(1),  np.maximum(((temp_out - scale_add) / scale_mult),1).astype(int) )

    # Convert nan back to 0
    temp_is_nan = is_nan
    if pad_rows_tmp>0: temp_is_nan = temp_is_nan[:-pad_rows_tmp]
    if pad_columns_tmp>0: temp_is_nan = temp_is_nan[:, :-pad_columns_tmp]
    
    temp_out[temp_is_nan] = 0

    # Make edits to the file
    band_ds.write(temp_out, 1)
    band_ds.close()
    
    del band_ds, image, image_AEC, image_R_surf, is_nan, R_conv, R_correction, R_correction_original_shape, temp_SR, temp_SR_water, temp_is_nan, temp_out
    gc.collect()
