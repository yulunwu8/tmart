# This file is part of TMart.
#
# Copyright 2024 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Overall control of AEC

def run(file, username, password, overwrite=False, AOT='MERRA2', n_photon=100_000, AOT_offset=0.0, n_jobs=100, mask_SWIR_threshold=None):
    '''Run adjacency-effect correction on satellite files. See 'Introduction - Adjacency-Effect Correction' for detailed instructions.
    
    Arguments:
        
    * ``file`` -- String. Path to satellite files. For L8/S2: provide path to the folder. For PRISMA: provide ACOLITE L1R file. SAFE.zip for S2 is supported: when doing so it is recommended to set ``overwrite`` as True to save space. 
    * ``username`` -- String. Username of EarthData account.
    * ``password`` -- String. Password of EarthData account.
    * ``overwrite`` -- Boolean. If overwrite the existing files. The default is False and it creates a folder in the same directory that starts with AEC in the name
    * ``AOT`` -- 'MERRA2': use ancillary data (default). Float: AOT at 550nm. 'NIR': calculate in T-Mart by finding dark pixels in NIR when considering the AE. 
    * ``n_photon`` -- Int. Number of photons in each T-Mart run, 100_000 is recommended for accurate results.
    * ``AOT_offset`` -- Float. Value added to AOT at 550nm. If resulted AOT is negative, it will be corrected to 0.
    * ``n_jobs`` -- Int. Number of jobs in Python multiprocessing. One CPU core processes one job at a time. n_photon is evenly distributed across the jobs.
    * ``mask_SWIR_threshold`` -- Float. Reflectance threshold in a SWIR band used to mask non-water pixels in the processing. If specified, this overwrites the value in config.txt. 

    Example usage::
        
        import tmart
        file = 'user/test/S2A_MSIL1C_20160812T143752_N0204_R096_T20MKB_20160812T143749.SAFE'
        username = 'abcdef'
        password = '123456'
        
        # T-Mart uses multiprocessing, which needs to be wrapped in 'if __name__ == "__main__":' for Windows users
        if __name__ == "__main__":
            tmart.AEC.run(file, username, password)

    '''
    
    import tmart
    import sys, os
    
    import time
    from importlib.metadata import version
    
    # Identify directories and files 
    file, file_is_dir = tmart.AEC.identify_input(file)
    home_folder = os.path.dirname(file) # where ancillary file and printed info are stored, default: the dir of the file
    basename = os.path.basename(file)
    basename_before_period = basename.split('.')[0]
    
    # If not overwrite, make a copy of the file(s) in the same directory  
    if not overwrite:
        basename_new = 'AEC_' + basename
        file_new = os.path.join(home_folder,basename_new)
        if os.path.exists(file_new):
            print('\nNew directory ' + str(file_new) + ' already exists, proceeds to AEC. \n')
        else: 
            # make the copy 
            print('\nCopying from ' + str(file) + ' to ' + str(file_new) + '... \n')
            if file_is_dir: 
                from shutil import copytree
                copytree(file, file_new)
            else:
                from shutil import copy
                copy(file, file_new)
        file = file_new
    
    # For S2/L89, ancillary file and printed info are stored in the image folder 
    if file_is_dir: home_folder = file

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
    print('AOT_offset: ' + str(AOT_offset))
    
    # S2/L89
    if file_is_dir: 
        tmart_out = tmart.AEC.run_regular(file, username, password, AOT, AOT_offset, n_photon, AEC_record, basename_before_period, n_jobs, mask_SWIR_threshold)
    # ACOLITE L1R 
    else:
        tmart_out = tmart.AEC.run_acoliteL1R(file, username, password, AOT, AOT_offset, n_photon, AEC_record, basename_before_period, n_jobs)
    
    # Stop logging 
    sys.stdout = orig_stdout

    return tmart_out
    