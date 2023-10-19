# This file is part of TMart.
#
# Copyright 2023 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Overall control of AEC

def run(file, username, password, overwrite=False, AOT = None, n_photon = 100_000):
    '''Run adjacency-effect correction on satellite files. Currently only supports Sentinel-2 MSI and Landsat 8 OLI products. See 'Introduction - Adjacency-Effect Correction' for details instructions.
    
    Arguments:
        
    * ``file`` -- String. Path to satellite files. 
    * ``username`` -- String. Username of EarthData account.
    * ``password`` -- String. Password of EarthData account.
    * ``overwrite`` -- Boolean. If overwrite the existing files. The default is False and it creates a folder in the same directory that starts with AEC in the name
    * ``AOT`` -- Float. AOT at 550nm, calculated in T-Mart if not specified.
    * ``n_photon`` -- Int. Number of photons in each T-Mart run, 100_000 is recommended for accurate results.   

    Example usage::
        
        file = 'user/test/S2A_MSIL1C_20160812T143752_N0204_R096_T20MKB_20160812T143749.SAFE'
        username = 'abcdef'
        password = '123456'
        tmart.AEC.run(file, username, password)

    '''
    
    import tmart
    import sys, os
    import time

    # If not overwrite, creates a new folder in the same directory  
    if not overwrite:
        from distutils.dir_util import copy_tree
        file_parts = file.rsplit('/',1)
        file_new = file_parts[0] + '/AEC_' + file_parts[1]
        if os.path.exists(file_new):
            print('\nNew directory ' + str(file_new) + ' already exists, proceeds to AEC. \n')
        else: 
            print('\nCopying from ' + str(file) + ' to ' + str(file_new) + '... \n')
            copy_tree(file, file_new)
        file = file_new

    # Check if AEC has been done, if so, exit 
    AEC_record = '{}/{}'.format(file,'AEC_completed.txt')
    if os.path.exists(AEC_record): 
        sys.exit('Record of adjacency correction found: ' + str(AEC_record) + 
                 '. Corrcting for the AE twice is not recommended and will lead to unpredictable results. ' + 
                 'If you’re confident that the imagery has not been corrected for the AE, delete this file to proceed.')

    # Start logging in txt file
    orig_stdout = sys.stdout
    log_file = 'tmart_log_' + time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime(time.time())) + '.txt'
    log_file = '{}/{}'.format(file,log_file)

    class Logger:
        def __init__(self, filename):
            self.console = sys.stdout
            self.file = open(filename, 'w')
            self.file.flush()
        def write(self, message):
            self.console.write(message)
            self.file.write(message)
        def flush(self):
            self.console.flush()
            self.file.flush()

    sys.stdout = Logger(log_file)
    print('System time: ' + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())))
    print('File: ' + str(file))
    
    # Read configuration
    config = tmart.AEC.read_config()

    # Identify sensor
    print('\nReading image files: ')
    print(file)
    sensor = tmart.AEC.identify_sensor(file)
    
    # Extract metadata
    if sensor == 'S2A' or sensor == 'S2B':
        metadata = tmart.AEC.read_metadata_S2(file,config, sensor)
    elif sensor == 'L8':
        metadata = tmart.AEC.read_metadata_L8(file,config)
    else: sys.exit('Warning: unrecognized sensor')
    metadata['sensor'] = sensor
    
    # Print metadata 
    print('Metadata: ')
    for k, v in metadata.items():
        print(str(k) + ': '  + str(v))
        
    # Get ancillary information from NASA Ocean Color
    anci = tmart.AEC.get_ancillary(metadata, username, password)
    
    # Compute cloud and non-Water masks 
    print('\nComputing masks: ')
    mask_cloud = tmart.AEC.compute_masks(metadata, config, 'cloud')
    mask_all   = tmart.AEC.compute_masks(metadata, config, 'all')   
    print('Done')
    
    # Estimate AOT using the 860nm band 
    if AOT == None:
        print('\nEstimating AOT from the NIR band: ')
        AOT = tmart.AEC.get_AOT(metadata, config, anci, mask_cloud, mask_all, n_photon)
    else:
        print('\nUser input AOT: ' + str(AOT))
    
    # Write atm information
    tmart.AEC.write_atm_info(file, anci, AOT)
    print('\nWrote aerosol and atmosphere information.')
    
    # Make a record file for AEC 
    file_AEC_record = open(AEC_record,"w")
    file_AEC_record.flush()
    
    # AEC for each of the specified bands 
    for i in range(len(metadata['AEC_bands_name'])):
        AEC_band_name = metadata['AEC_bands_name'][i] 
        AEC_band_6S = metadata['AEC_bands_6S'][i]
        wl = metadata['AEC_bands_wl'][i]
        print('\n============= AEC: {} ==================='.format(AEC_band_name))
        tmart.AEC.AEC(AEC_band_name, AEC_band_6S, wl, AOT, metadata, config, anci, mask_cloud, mask_all, n_photon)
        file_AEC_record.write(str(AEC_band_name) + '\n')
        file_AEC_record.flush()
        
    # Close AEC record 
    file_AEC_record.close()
    
    # Stop logging 
    sys.stdout = orig_stdout

    return 0
    
    