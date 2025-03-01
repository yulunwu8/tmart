# This file is part of TMart.
#
# Copyright 2024 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Run on Landsat and Sentinel series 

def run_regular(file, username, password, AOT, AOT_offset, n_photon, AEC_record, basename, njobs, mask_SWIR_threshold): 
 
    import tmart
    import sys, os
    
    # read configuration
    config = tmart.AEC.read_config(mask_SWIR_threshold)

    # identify sensor
    print('\nReading image files: ')
    print(file)
    sensor = tmart.AEC.identify_sensor(file)
    
    # extract metadata
    if sensor == 'S2A' or sensor == 'S2B' or sensor == 'S2C':
        metadata = tmart.AEC.read_metadata_S2(file, config, sensor)
    elif sensor == 'L8' or sensor == 'L9':
        metadata = tmart.AEC.read_metadata_Landsat(file, config, sensor)
    else: sys.exit('Warning: unrecognized sensor')
    
    metadata['sensor'] = sensor
    
    # print metadata 
    print('Metadata: ')
    for k, v in metadata.items():
        print(str(k) + ': '  + str(v))
        
    # get ancillary information from NASA Ocean Color
    anci = tmart.AEC.get_ancillary(metadata, username, password)
    
    # compute cloud and non-Water masks 
    print('\nComputing masks: ')
    mask_cloud = tmart.AEC.compute_masks(metadata, config, 'cloud')
    mask_all   = tmart.AEC.compute_masks(metadata, config, 'all')   
    print('Done')
    
    # plot water-extent preview
    if not config['AE_land']=='True':
        print('\nPlotting water-extent preview...')
        tmart.AEC.plot_water_extent(metadata, config, mask_all)
        
    # AOT
    if AOT == 'NIR':
        print('\nEstimating AOT from the NIR band: ')
        AOT = tmart.AEC.get_AOT(metadata, config, anci, mask_cloud, mask_all, n_photon)
    elif AOT == 'MERRA2':
        AOT = anci['AOT_MERRA2']
        print('\nUsing AOT550 from MERRA2: ' + str(AOT))
    else:
        print('\nUser input AOT: ' + str(AOT))
    
    # add offset
    AOT = max(0, AOT + AOT_offset)
    
    # write atm information
    tmart.AEC.write_atm_info(file, basename, anci, AOT)
    print('\nWrote aerosol and atmosphere information.')
    
    # make a record file for AEC 
    file_AEC_record = open(AEC_record,"w")
    file_AEC_record.flush()
    
    # AEC for each of the specified bands 
    for i in range(len(metadata['AEC_bands_name'])):
        AEC_band_name = metadata['AEC_bands_name'][i] 
        AEC_band_6S = metadata['AEC_bands_6S'][i]
        wl = metadata['AEC_bands_wl'][i]
        print('\n============= AEC: {} ==================='.format(AEC_band_name))
        tmart.AEC.AEC(AEC_band_name, AEC_band_6S, wl, AOT, metadata, config, anci, mask_cloud, mask_all, n_photon, njobs)
        file_AEC_record.write(str(AEC_band_name) + '\n')
        file_AEC_record.flush()
        
    # close AEC record 
    file_AEC_record.close()
    
    return 0
    