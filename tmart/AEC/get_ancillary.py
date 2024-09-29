# This file is part of TMart.
#
# Copyright 2024 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# Extract aerosol type, ozone and water vapour from ancillary data from NASA Ocean Color 

def get_ancillary(metadata, username, password):
    
    import tmart
    
    # List all files to download 
    files = tmart.AEC.anci_list_files(metadata)
    
    # Download the files  
    file = metadata['file']
    
    print('\nDownloading ancillary aerosol data from NASA Ocean Color: ')
    file_aer_h0 = tmart.AEC.anci_download(file, files['AER'][0], username, password)
    file_aer_h1 = tmart.AEC.anci_download(file, files['AER'][1], username, password)
    
    print('\nDownloading ancillary meteorological data from NASA Ocean Color: ')
    file_met_h0 = tmart.AEC.anci_download(file, files['MET'][0], username, password)
    file_met_h1 = tmart.AEC.anci_download(file, files['MET'][1], username, password)
    
    # Extract aerosol type
    r_maritime = tmart.AEC.anci_get_AER([file_aer_h0,file_aer_h1],metadata)
    
    # Extract ozone and water vapour 
    OWV = tmart.AEC.anci_get_OWV([file_met_h0,file_met_h1],metadata)
    
    return {'r_maritime': r_maritime[0],
            'Angstrom_exp': r_maritime[1],
            'SSA': r_maritime[2],
            'AOT_MERRA2': r_maritime[3],
            'ozone': OWV[0],
            'water_vapour': OWV[1]}
    
