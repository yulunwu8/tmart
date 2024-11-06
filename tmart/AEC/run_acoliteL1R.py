# This file is part of TMart.
#
# Copyright 2024 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Run on ACOLITE L1R files, currently supports PRSIMA only 

def run_acoliteL1R(file, username, password, AOT, AOT_offset, n_photon, AEC_record, basename, njobs):
 
    import tmart
    import sys, os, time
    import netCDF4 as nc4
    import numpy as np
    from scipy import signal, interpolate
        
    # Read configuration
    config = tmart.AEC.read_config()
    
    # Open netcdf file 
    dset = nc4.Dataset(file, 'r+')
    
    # wavelengths and widths 
    WLs = dset.getncattr('band_waves')
    WWs = dset.getncattr('band_widths')
    
    ### Set up metadata  
    window_size = 101
    reshape_factor = 5

    metadata = {}
    
    # Dimension 
    metadata['height'] = 1000
    metadata['width'] = 1000

    # Reshaped dimension for AEC
    height_reshaped = int(metadata['height'] / reshape_factor)
    width_reshaped = int(metadata['width'] / reshape_factor)
    
    # metadata and time 
    metadata['file'] = os.path.dirname(file)
    metadata['time'] = dset.getncattr('isodate')[0:19]
    
    metadata['sza'] = dset.getncattr('sza')
    metadata['saa'] = dset.getncattr('saa')
    metadata['vza'] = 0 # dset.getncattr('vza') # missing vaa for now so we assume 0 
    metadata['vaa'] = 0 # dset.getncattr('vaa')
    
    # find mean coordinates and direction to north 
    metadata['lat'] = np.mean(dset['lat'][:])
    metadata['lon'] = np.mean(dset['lon'][:])
    
    # Up relative to north 
    tm_PRISMA_up = tmart.AEC.read_PRISMA_north(file)
    sun_dir = [metadata['sza'], (metadata['saa']+270 - tm_PRISMA_up)%360]
        
    # Print metadata 
    print('Metadata: ')
    for k, v in metadata.items():
        print(str(k) + ': '  + str(v))
    
    # Get ancillary information from NASA Ocean Color
    anci = tmart.AEC.get_ancillary(metadata, username, password)
    
    
    ### Masks 
    
    # Cloud mask, may need to find band closest to 1373 in case of future wavelength shift 
    image_cloud = dset['rhot_1373'][:]
    image_cloud_compressed = image_cloud.reshape([height_reshaped, reshape_factor, width_reshaped, reshape_factor]).mean(3).mean(1)
    mask_cloud_compressed = image_cloud_compressed > float(config['mask_cirrus_threshold'])
    
    # SWIR, may need to find band closest to 1600 in case of future wavelength shift 
    mask_SWIR = dset['rhot_1596'][:] > float(config['mask_SWIR_threshold'])
    
    # TOA out of limit, all bands > 0.3 
    mask_highTOA = np.zeros(np.shape(mask_SWIR))    
    for WL in WLs: 
        print('High TOA: {:.1f} nm'.format(WL))
        rhot_wl = 'rhot_' + str(int(round(WL)))
        image_highTOA = dset[rhot_wl][:]
        if np.ma.count_masked(image_highTOA) > 0: 
            sys.exit('T-Mart assumed PRISMA has no masked values, please use another image or contact Yulun at yulunwu8@gmail.com')
        mask_highTOA = np.logical_or(mask_highTOA, image_highTOA > float(config['mask_highTOA_threshold']))
    
    # Combine highTOA and SWIR masks, original resolution 
    mask_all = np.logical_or(mask_highTOA, mask_SWIR)
    
    
    # AOT
    if AOT == 'NIR':
        sys.exit('The default estimating AOT from the NIR band is not tested on PRISMA data, please use MERRA2 or user input instead.')
    elif AOT == 'MERRA2':
        AOT = anci['AOT_MERRA2']
        print('\nUsing AOT from MERRA2: ' + str(AOT))
    else:
        print('\nUser input AOT: ' + str(AOT))
    
    # Add offset
    AOT = max(0, AOT + AOT_offset)
    
    # Write atm information
    tmart.AEC.write_atm_info(metadata['file'], basename, anci, AOT)
    print('\nWrote aerosol and atmosphere information.')
    
    
    ### T-Mart calculate correction parameters 

    # Loop through 5 wavelengths, generating 4 look-up tables    
    tm_wls_interp = np.array([406.9,550,783,1000,2429])
    list_conv_window_1  = []
    list_F_correction   = []
    list_F_captured     = []
    list_R_atm          = []
        
    for tm_wl_interp in tm_wls_interp: 
        wl = tm_wl_interp
        AEC_parameters = tmart.AEC.get_parameters(n_photon = n_photon, wl = wl, 
                                                  target_pt_direction=[180,0], sun_dir=sun_dir, 
                                                  atm_profile = anci, 
                                                  aerosol_type = anci['r_maritime'], aot550 = AOT, 
                                                  cell_size = 30*reshape_factor,
                                                  window_size = window_size, njobs=njobs)
        
        conv_window_1   = AEC_parameters['conv_window_1']
        F_correction    = AEC_parameters['F_correction']
        F_captured      = AEC_parameters['F_captured']
        R_atm           = AEC_parameters['R_atm']
    
        list_conv_window_1.append(conv_window_1)
        list_F_correction.append(F_correction)
        list_F_captured.append(F_captured)
        list_R_atm.append(R_atm)
    
    # make np arrays 
    list_conv_window_1 = np.array(list_conv_window_1)
    list_F_correction = np.array(list_F_correction)
    list_F_captured = np.array(list_F_captured)
    list_R_atm = np.array(list_R_atm)
    
    # interpolation objects  
    interp_conv_window_1 = interpolate.interp1d(tm_wls_interp, list_conv_window_1, kind='linear', axis=0)
    interp_F_correction = interpolate.interp1d(tm_wls_interp, list_F_correction, kind='linear', axis=0)
    interp_F_captured = interpolate.interp1d(tm_wls_interp, list_F_captured, kind='linear', axis=0)
    interp_R_atm = interpolate.interp1d(tm_wls_interp, list_R_atm, kind='linear', axis=0)
        
    # Make a record file for AEC 
    file_AEC_record = open(AEC_record,"w")
    file_AEC_record.flush()
    
    # Loop through all bands 
    for i in range(len(WLs)):
        
        # Wavelength and FWHM
        WL = WLs[i]
        WW = WWs[i]
        
        print('\n----------')
        print('Processing wavelength {:.2f} nm'.format(WL))
        
        # If tgas < 0.7, then skip AEC for this wavelength 
        tgas = tmart.AEC.compute_gas_transmittance(metadata, anci, WL, WW)
        if tgas < 0.7: 
            print('Skipping, tgas: {:.2f} < 0.7'.format(tgas))
            continue
        else:
            print('tgas: {:.2f}'.format(tgas))
            
        # Skip if out of range 
        if WL > tm_wls_interp.max():
            print('Skipping band out of correction-parameter range')
            continue
    
        # Variable name 
        rhot_wl = 'rhot_' + str(int(round(WL)))
        
        # interpolate AEC parameters 
        conv_window_1   = interp_conv_window_1(WL)
        F_correction    = interp_F_correction(WL)
        F_captured      = interp_F_captured(WL)
        R_atm           = interp_R_atm(wl)
        
        # Read image 
        image = dset[rhot_wl][:].data
        
        # Turn negative TOA to 0 
        image[image<0] = 0
        
        # Reshape
        image_AEC = image.reshape([height_reshaped, reshape_factor, 
                                   width_reshaped, reshape_factor]).mean(3).mean(1)
        
        # image_R_surf: reshaped 
        image_R_surf = image_AEC - R_atm
        print('\nNumber of negative R_surf pixels: ' + str(np.sum(image_R_surf < 0)) + '/' + str(height_reshaped * width_reshaped))
        
        # Turn negative to 0
        image_R_surf[image_R_surf<0] = 0
        
        # Read AEC parameters 
        hp_cloud_contribution = float(config['cloud_contribution'])
        hp_AE_land = config['AE_land'] == 'True'
        
        # Reduce contribution of clouds according to setting 
        image_R_surf[mask_cloud_compressed] = image_R_surf[mask_cloud_compressed] * hp_cloud_contribution
        
        # Convolution 
        start_time = time.time()
        print("\nConvolution started ")
        filter_kernel = np.flip(conv_window_1) # it's flipped in convolve by default
        R_conv = signal.convolve2d(image_R_surf, filter_kernel, mode='same', boundary='fill', fillvalue=image_R_surf.mean())
        print("Convolution completed: %s seconds " % (time.time() - start_time))
        
        
        # Smoothing the edges 
        if reshape_factor >1 and not hp_AE_land:
            image_water = np.where(mask_all, np.nan, image)
            
            # Suppress warning of mean of empty slice 
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                image_water_AEC = image_water.reshape([height_reshaped, reshape_factor, width_reshaped, reshape_factor])
                image_water_AEC = np.nanmean(image_water_AEC,3)
                image_water_AEC = np.nanmean(image_water_AEC,1)
        
            # with non-water as nan 
            image_water_R_surf = image_water_AEC - R_atm
            
            # Edit image_R_surf: where there's water image_R_surf, keep it, where there's nan, use original image_R_surf
            # image_R_surf can only be edited after convolution 
            image_R_surf = np.where(np.isnan(image_water_R_surf), image_R_surf, image_water_R_surf)   
            
        # Correction to make 
        R_correction = (R_conv - image_R_surf) * F_correction 
        print('\nNumber of pixels where R_correction > R_surf : ' + str(np.sum(R_correction>image_R_surf)) + '/' + str(height_reshaped * width_reshaped))
    
        # Back to the original size 
        R_correction_original_shape = np.repeat(np.repeat(R_correction, reshape_factor, axis=0), reshape_factor, axis=1)
    
        # Smoothing the gridline artifacts 
        if reshape_factor>1: 
        
            # Keep only water, short name for R_correction_original_shape_water
            RC_water = np.where(mask_all, np.nan, R_correction_original_shape)
            RC_water_fill_nan = tmart.AEC.fillnan(RC_water)
            
            filter_size = 3
            filter_kernel = np.full((filter_size,filter_size), 1/(filter_size**2))
            
            RC_water_smooth = signal.convolve2d(RC_water_fill_nan, filter_kernel, mode='same', boundary='fill', fillvalue=np.nanmean(RC_water))
            
            # place RC_water_smooth in R_correction_original_shape where there's water
            R_correction_original_shape = np.where(mask_all, R_correction_original_shape, RC_water_smooth)
    
        # remove AE from TOA reflectance 
        temp_SR = image - R_correction_original_shape - R_atm
        
        # Correcting for non-linearity of the ratio of environmental irradiance to surface reflectance for homogeneous Lambertian surfaces
        # This is implemented for water only
        # Input has to include land, because the function uses the env irradiance of the average reflectance across the scene
        print("\nIrradiance correction: ")
        temp_SR_water = tmart.AEC.irradiance_correction(image = temp_SR, wl_RC = wl/1000, band = None,
                                                        tm_vza = metadata['vza'], tm_vaa = metadata['vaa'], 
                                                        tm_sza = metadata['sza'], tm_saa = metadata['saa'],
                                                        atm_profile = anci, 
                                                        aerosol_type = anci['r_maritime'], aot550 = AOT)
        
        # If land correction is needed 
        if hp_AE_land:
            temp_out = np.where(mask_all, temp_SR, temp_SR_water)
        else:
            temp_out = np.where(mask_all, (image-R_atm), temp_SR_water)

        temp_out = temp_out + R_atm
        
        # Modify L1R NC file 
        dset[rhot_wl][:] = temp_out
            
        # Write AEC record
        file_AEC_record.write(str(rhot_wl) + '\n')
        file_AEC_record.flush()
    
    # Close AEC record 
    file_AEC_record.close()
    
    # Close file 
    dset.close()
    
    return 0
    
 