# This file is part of TMart.
#
# Copyright 2023 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# Read S2 metadata

def read_metadata_S2(file,config,sensor):
    import tmart
    import os, sys
    import mgrs
    import Py6S, math
    
    files = os.listdir(file)
    metadata = {}
    metadata['file'] = file
    
    # resolution for masks 
    mask_res = [10,20,60] # mask resolution 
    mask_res.append(int( 10 *  int(config['reshape_factor_S2'])))
    metadata['mask_res'] = mask_res
    
    # high TOA bands for masks 
    highTOA_band_names = ['B01','B02','B03','B04','B05','B06','B07','B08','B8A','B11','B12']
    highTOA_band_names.remove(config['S2_SWIR_band'])
    metadata['highTOA_band_names'] = highTOA_band_names
    
    # masks 
    metadata['cirrus_mask'] = config['S2_cirrus_band']
    metadata['SWIR_mask'] = config['S2_SWIR_band']
    
    # Bands to be AECed
    metadata['AEC_bands_name'] = ['B01','B02','B03','B04','B05','B06','B07','B08','B8A','B11','B12']
    
    if sensor == 'S2A': 
        metadata['AEC_bands_6S'] = [Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_01),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_02),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_03),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_04),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_05),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_06),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_07),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_08),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_8A),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_11),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_12)]
        metadata['AEC_bands_wl'] = [442.7, 492.4, 559.8, 664.6, 704.1, 740.5, 782.8, 832.8, 864.7, 1613.7, 2202.4]
    elif sensor == 'S2B':
        metadata['AEC_bands_6S'] = [Py6S.Wavelength(Py6S.PredefinedWavelengths.S2B_MSI_01),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2B_MSI_02),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2B_MSI_03),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2B_MSI_04),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2B_MSI_05),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2B_MSI_06),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2B_MSI_07),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2B_MSI_08),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2B_MSI_8A),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2B_MSI_11),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2B_MSI_12)]
        metadata['AEC_bands_wl'] = [442.2, 492.1, 559.0, 664.9, 703.8, 739.1, 779.7, 832.9, 864.0, 1610.4, 2185.5]        
        
    # Find central coordinates of MGRS tile 
    m = mgrs.MGRS()
    
    for i, fname in enumerate(files):
        tmp = fname.split('.')
        path = '{}/{}'.format(file,fname)
        
        # Granules
        if (fname == 'GRANULE'):
            granules = os.listdir(path)
            
            # Check if there is only one granule file 
            n_granule = 0
            
            for granule in granules:
                if granule[0]=='.':continue
                
                n_granule += 1
                if n_granule>1: sys.exit('Warning: more than 1 granule')
                
                metadata['granule'] = '{}/{}/{}/IMG_DATA/'.format(file,fname,granule)
                metadata['MGRS_tile'] = granule.split('_')[1][1:]
                metadata['QI_DATA'] = '{}/{}/{}/QI_DATA'.format(file,fname,granule)
                
                # MGRS
                tile = metadata['MGRS_tile'] + '54905490'
                d = m.toLatLon(tile)
                metadata['lat'] = d[0]
                metadata['lon'] = d[1]
                
                # Band files 
                image_files = os.listdir(metadata['granule'])
                for image in image_files: 
                    if image[0]=='.':continue
                    if image[-4:]=='.xml':continue
                    tmp = image.split('_')
                    metadata[tmp[-1][0:3]] = '{}/{}/{}/IMG_DATA/{}'.format(file,fname,granule,image)
                
                # Read scene metadata  
                path = '{}/{}/{}/'.format(file,fname,granule)
                granule_files = os.listdir(path)
                for j, grfname in enumerate(granule_files):
                    tmp = grfname.split('.')
                    path = '{}/{}/{}/{}'.format(file,fname,granule,grfname)
                    if (tmp[-1] == ('xml')) & ('MTD' in tmp[0]):
                        xml = tmart.AEC.read_xml_S2(path)
                        metadata.update(xml)
    
    # Scaling factors 
    xml2 = tmart.AEC.read_xml_S2_scene('{}/MTD_MSIL1C.xml'.format(file))
    metadata.update(xml2)
    
    # Others 
    metadata['resolution'] = 10
    metadata['reshape_factor'] = int(config['reshape_factor_S2'])
    if metadata['reshape_factor'] % 6 > 0: sys.exit('Warning: reshape_factor_S2 in config must be divisible by 6.')
    
    metadata['window_size'] = int(config['window_size'])
    metadata['AEC_resolution'] = int(config['reshape_factor_S2']) * 10 # resolution 

    metadata['height'] = 10_980
    metadata['width'] = 10_980
    
    # The smallest number divisible by both 6 and reshape_factor
    # 6 comes from 60m/10m resolution, remainder of 6*reshape_f divided by greatest common divisor between 6 and reshape_f
    temp_reshape_factor = (metadata['reshape_factor'] * 6) // math.gcd(metadata['reshape_factor'], 6)
    
    # Height and width for AEC, with a few extra rows and columns 
    metadata['AEC_height'] = math.ceil(metadata['height'] / temp_reshape_factor) * temp_reshape_factor
    metadata['AEC_width'] = math.ceil(metadata['width'] / temp_reshape_factor) * temp_reshape_factor
    
    return metadata 
    
    
    
    