# This file is part of TMart.
#
# Copyright 2023 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# Read L8/9 metadata

def read_metadata_L8(file, config):
    import pandas as pd
    import glob, os, Py6S
    import math 
    
    # MTL file 
    mtl = glob.glob('{}/{}'.format(file, '*MTL.txt'))[0]
    df = pd.read_csv(mtl, delimiter = '=', skipinitialspace=True)
    
    # Get the first and second columns by integer positions
    first_column = df.iloc[:, 0]
    first_column = first_column.str.replace(' ', '')
    second_column = df.iloc[:, 1]

    # Convert to dictionary
    result_dict = dict(zip(first_column, second_column))

    # Extract values 
    SUN_AZIMUTH = float(result_dict['SUN_AZIMUTH'])
    SUN_ELEVATION = float(result_dict['SUN_ELEVATION'])
    time = result_dict['DATE_ACQUIRED']+'T'+result_dict['SCENE_CENTER_TIME']
    lat = (float(result_dict['CORNER_UL_LAT_PRODUCT']) + float(result_dict['CORNER_UR_LAT_PRODUCT']) + 
           float(result_dict['CORNER_LL_LAT_PRODUCT']) + float(result_dict['CORNER_LR_LAT_PRODUCT']))/4
    lon = (float(result_dict['CORNER_UL_LON_PRODUCT']) + float(result_dict['CORNER_UR_LON_PRODUCT']) + 
           float(result_dict['CORNER_LL_LON_PRODUCT']) + float(result_dict['CORNER_LR_LON_PRODUCT']))/4
    tm_sun_dir=[90.0-SUN_ELEVATION, (SUN_AZIMUTH+270)%360] # mean_saa = 0 => sun is in the north => 270 in T-Mart
    
    # Create dictionary 
    metadata =  {'file': file,
                 'granule': file,
                 'lat': lat,
                 'lon': lon,
                 'time':time,
                 'saa':SUN_AZIMUTH,
                 'sza':90.0-SUN_ELEVATION,
                 'vaa':0,
                 'vza':0,
                 'tm_pt_dir': [180,0],
                 'tm_sun_dir': tm_sun_dir}
    
    # Scaling factors 
    metadata['B1_mult'] = float(result_dict['REFLECTANCE_MULT_BAND_1']) 
    metadata['B2_mult'] = float(result_dict['REFLECTANCE_MULT_BAND_2']) 
    metadata['B3_mult'] = float(result_dict['REFLECTANCE_MULT_BAND_3']) 
    metadata['B4_mult'] = float(result_dict['REFLECTANCE_MULT_BAND_4']) 
    metadata['B5_mult'] = float(result_dict['REFLECTANCE_MULT_BAND_5']) 
    metadata['B6_mult'] = float(result_dict['REFLECTANCE_MULT_BAND_6']) 
    metadata['B7_mult'] = float(result_dict['REFLECTANCE_MULT_BAND_7']) 
    metadata['B9_mult'] = float(result_dict['REFLECTANCE_MULT_BAND_9']) 
    metadata['B1_add'] = float(result_dict['REFLECTANCE_ADD_BAND_1']) 
    metadata['B2_add'] = float(result_dict['REFLECTANCE_ADD_BAND_2']) 
    metadata['B3_add'] = float(result_dict['REFLECTANCE_ADD_BAND_3']) 
    metadata['B4_add'] = float(result_dict['REFLECTANCE_ADD_BAND_4']) 
    metadata['B5_add'] = float(result_dict['REFLECTANCE_ADD_BAND_5']) 
    metadata['B6_add'] = float(result_dict['REFLECTANCE_ADD_BAND_6']) 
    metadata['B7_add'] = float(result_dict['REFLECTANCE_ADD_BAND_7'])
    metadata['B9_add'] = float(result_dict['REFLECTANCE_ADD_BAND_9'])
    
    # Band and mask files 
    files = os.listdir(file)
    for image in files:
        if image[0]=='.':continue
        tmp = image.split('_')
        if tmp[-1][0] == 'B' and tmp[-1][-3:] == 'TIF':
            metadata[tmp[-1].split('.')[0]]  = os.path.join(file,image)
        elif tmp[-1] == 'PIXEL.TIF': 
            metadata['cloud_mask'] = os.path.join(file,image)
        else:
            pass
    
    # Resolution for masks 
    metadata['mask_res'] = [30, int( 30 *  int(config['reshape_factor_L8']))]
    
    # High TOA bands for masks 
    highTOA_band_names = ['B1','B2','B3','B4','B5','B6','B7'] # bands to use in the highTOA mask 
    highTOA_band_names.remove(config['L8_SWIR_band'])
    metadata['highTOA_band_names'] = highTOA_band_names

    # Masks 
    metadata['cirrus_mask'] = config['L8_cirrus_band']
    metadata['SWIR_mask'] = config['L8_SWIR_band']
  
    # Bands to be AECed
    metadata['AEC_bands_name'] = ['B1','B2','B3','B4','B5','B6','B7']
    metadata['AEC_bands_6S'] = [Py6S.Wavelength(Py6S.PredefinedWavelengths.LANDSAT_OLI_B1),
                                Py6S.Wavelength(Py6S.PredefinedWavelengths.LANDSAT_OLI_B2),
                                Py6S.Wavelength(Py6S.PredefinedWavelengths.LANDSAT_OLI_B3),
                                Py6S.Wavelength(Py6S.PredefinedWavelengths.LANDSAT_OLI_B4),
                                Py6S.Wavelength(Py6S.PredefinedWavelengths.LANDSAT_OLI_B5),
                                Py6S.Wavelength(Py6S.PredefinedWavelengths.LANDSAT_OLI_B6),
                                Py6S.Wavelength(Py6S.PredefinedWavelengths.LANDSAT_OLI_B7)]
    
    # Center wavelength, source: http://gsics.atmos.umd.edu/pub/Development/20171106/5k_Ong_Landsat8_Lunar_Calibrations.pdf
    metadata['AEC_bands_wl'] = [443, 482, 561.4, 654.6, 864.7, 1608.9, 2200.7]
    
    # Others 
    metadata['resolution'] = 30
    metadata['reshape_factor'] = int(config['reshape_factor_L8'])
    metadata['window_size'] = int(config['window_size'])
    metadata['AEC_resolution'] = int(config['reshape_factor_L8']) * 30 # resolution 
    metadata['height'] = int(result_dict['REFLECTIVE_LINES']) 
    metadata['width'] = int(result_dict['REFLECTIVE_SAMPLES']) 
    
    # Height and width for AEC, with a few extra rows and columns 
    metadata['AEC_height'] = math.ceil(metadata['height'] / int(config['reshape_factor_L8'])) * int(config['reshape_factor_L8'])
    metadata['AEC_width'] = math.ceil(metadata['width'] / int(config['reshape_factor_L8'])) * int(config['reshape_factor_L8'])

    return metadata
    
 