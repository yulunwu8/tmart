# This file is part of TMart.
#
# Copyright 2023 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Overall control of AEC

def run(file, username, password, overwrite=False, AOT = 'MERRA2', n_photon = 100_000):
    '''Run adjacency-effect correction on satellite files. Currently only supports Sentinel-2 MSI and Landsat 8 OLI products. See 'Introduction - Adjacency-Effect Correction' for detailed instructions.
    
    Arguments:
        
    * ``file`` -- String. Path to satellite files. L8/S2: provide path to the folder. PRISMA: ACOLITE L1R file. 
    * ``username`` -- String. Username of EarthData account.
    * ``password`` -- String. Password of EarthData account.
    * ``overwrite`` -- Boolean. If overwrite the existing files. The default is False and it creates a folder in the same directory that starts with AEC in the name
    * ``AOT`` -- Float. AOT at 550nm, calculated in T-Mart if not specified. 'MERRA2' ancillary data is another option. 
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
    from importlib.metadata import version
    
    
    # Identify if directory or file, for S2/L8
    if os.path.isdir(file):
        file_is_dir = True
        home_folder = file
    # For ACOLITE L1R files, currently supports PRISMA data only 
    elif os.path.isfile(file) and file[-6:]=='L1R.nc': 
        file_is_dir = False
        home_folder = os.path.dirname(file)
    else: 
        sys.exit('File or folder cannot be identified: {}'.format(file))
    
    basename = os.path.basename(file)
    basename_before_period = basename.split('.')[0]
    
    # If not overwrite, make a copy of the file(s) in the same directory  
    if not overwrite:
        basename_new = 'AEC_' + basename
        file_new = os.path.join(home_folder,basename_new)
        if os.path.exists(file_new):
            print('\nNew directory ' + str(file_new) + ' already exists, proceeds to AEC. \n')
        else: 
            print('\nCopying from ' + str(file) + ' to ' + str(file_new) + '... \n')
            if file_is_dir: 
                from distutils.dir_util import copy_tree
                copy_tree(file, file_new)
            else:
                from shutil import copy
                copy(file, file_new)
        file = file_new

    # Check if AEC has been done, if so, exit 
    AEC_record = '{}/AEC_completed_{}.txt'.format(home_folder,basename_before_period)
    if os.path.exists(AEC_record): 
        sys.exit('Record of adjacency correction found: ' + str(AEC_record) + 
                 '. Corrcting for the AE twice is not recommended and will lead to unpredictable results. ' + 
                 'If you’re confident that the imagery has not been corrected for the AE, delete this file to proceed.')

    # Start logging in txt file
    orig_stdout = sys.stdout
    log_file = 'tmart_log_' + time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime(time.time())) + '.txt'
    log_file = '{}/{}'.format(home_folder,log_file)

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
    print('T-Mart version: ' + str(version('tmart')))
    print('System time: ' + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time())))
    print('File: ' + str(file))
    
    # S2/L8
    if file_is_dir: 
        tmart_out = tmart.AEC.run_regular(file, username, password, AOT, n_photon, AEC_record, basename_before_period)
    # ACOLITE L1R 
    else:
        tmart_out = tmart.AEC.run_acoliteL1R(file, username, password, AOT, n_photon, AEC_record, basename_before_period)
    
    # Stop logging 
    sys.stdout = orig_stdout

    return tmart_out
    
    