# This file is part of TMart.
#
# Copyright 2024 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# Create a masks at different resolutions 

def compute_masks(metadata, config, mask_type):
    # mask_type: 'cloud' or 'all' 
    # output masks in a dictionary, with resolution as keys. E.g., masks['10m'], masks['20m'], etc.
    
    import rasterio, sys, os, gc
    import numpy as np
    
    sensor = metadata['sensor']
    
    # Input from metadata 
    mask_res = metadata['mask_res']
    highTOA_band_names = metadata['highTOA_band_names']
    
    # Number of columns and rows to pad
    pad_columns = metadata['AEC_width'] - metadata['width'] 
    pad_rows = metadata['AEC_height'] - metadata['height'] 
    
    # Mask function 
    def mask_threshold(band_name, threshold, mask_NAN, reshape = True): 
        # Output either an array or a dictionary, mask['10m'], etc. 
        
        # Read file and array 
        band_file = metadata[band_name]
        band_ds = rasterio.open(band_file)
        band_array = band_ds.read(1)
        
        # Resolution of this particular band 
        res_band = int(abs(band_ds.transform[0]))
        
        # Pad less for lower resolution 
        pad_rows_tmp = int(pad_rows/(res_band/metadata['resolution']))
        pad_columns_tmp = int(pad_columns/(res_band/metadata['resolution']))
        
        # Padding 
        band_array = np.pad(band_array, ((0, pad_rows_tmp), (0, pad_columns_tmp)), mode='constant', constant_values=0)
        
        # Scaling  
        scale_mult = metadata[str(band_name) + '_mult']
        scale_add = metadata[str(band_name) + '_add']
        mask_threshold = (float(config[threshold])  - scale_add ) / scale_mult
        
        # Make mask. If mask_NAN, then mask the NAN values 
        if mask_NAN: mask = np.logical_or(band_array > mask_threshold, band_array == 0) # In band_array,0 is masks 
        else: mask = band_array > mask_threshold
        
        # Add other resolution 
        if reshape:
            
            # Mask 0 as nan 
            band_array = band_array * 1.0
            is_nan = band_array == 0
            band_array[is_nan] = np.nan
            
            masks = {}
            
            # For loop of target resolution  
            for res in mask_res:

                # Height and width at this resolution 
                height_reshaped = int( metadata['AEC_height'] / (res / metadata['resolution'])) 
                width_reshaped = int( metadata['AEC_width'] / (res / metadata['resolution'])) 
    
                # If this is the resolution 
                if res == res_band:
                    masks[str(res) + 'm'] = mask
                    
                # If target resolution is coarser than band resolution, we take the reshaped means 
                elif res > res_band:
    
                    # Suppress warning of mean of empty slice 
                    import warnings
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", category=RuntimeWarning)
                        band_array_reshaped = band_array.reshape([height_reshaped, int(res / res_band), 
                                                                  width_reshaped, int(res / res_band)])
                        band_array_reshaped = np.nanmean(band_array_reshaped,3)
                        band_array_reshaped = np.nanmean(band_array_reshaped,1)
            
                    if mask_NAN: mask_reshaped = np.logical_or(band_array_reshaped > mask_threshold, np.isnan(band_array_reshaped))
                    else: mask_reshaped = band_array_reshaped > mask_threshold
                        
                    masks[str(res) + 'm'] = mask_reshaped
                
                # If target resolution is finer than band resolution, we repeat cells 
                elif res < res_band: 
                    
                    # stretch to the needed shape, 10980 by 10980 for S2 
                    mask_temp = np.repeat(np.repeat(mask, int(res_band/res), axis=0), int(res_band/res), axis=1)
                    masks[str(res) + 'm'] = mask_temp
                
            return masks
        
        else:
            return mask
    
    # Cirrus band 
    mask_cirrus = mask_threshold(metadata['cirrus_mask'], threshold = 'mask_cirrus_threshold', mask_NAN = False)
    
    # Cloud 
    if mask_type == 'cloud':
        print('Computing cloud mask...')
        
        # Sentinel-2
        if sensor == 'S2A' or sensor == 'S2B' or sensor == 'S2C':
            gml_file = "{}/MSK_CLOUDS_B00.gml".format(metadata['QI_DATA'])
            jp2_file = "{}/MSK_CLASSI_B00.jp2".format(metadata['QI_DATA'])
            
            # For imagery before processing baseline 4: Jan 25, 2022
            if os.path.exists(gml_file): 
        
                # Built-in cloud mask 
                import geopandas as gpd
                from rasterio.features import geometry_mask
                
                # Load a raster as the base of the mask 
                image = metadata['B02']
                
                with rasterio.open(image) as src:
                    # Read the raster data and transform
                    raster_data = src.read(1)
                    transform = src.transform
                    crs = src.crs
                
                try: 
                    # Read GML file 
                    gdf = gpd.read_file(gml_file)
            
                    # Create a mask using the GML polygons and the GeoTIFF metadata
                    mask_cloud = geometry_mask(gdf['geometry'], transform=transform, out_shape=raster_data.shape, invert=True)
                
                # Sometimes the GML file contains no information, assume no clouds in such case
                except:
                    mask_cloud = np.zeros_like(raster_data)
                
            # For imagery processing baseline 4
            elif os.path.exists(jp2_file): 
                band_ds = rasterio.open(jp2_file)
                band_array = band_ds.read(1)
                mask_cloud = band_array == 1
                mask_cloud = np.repeat(np.repeat(mask_cloud, 6, axis=0), 6, axis=1)
                
            else:
                sys.exit('Warning: cloud mask missing in {}.'.format(metadata['QI_DATA']))
                
        # Landsat 8
        elif sensor == 'L8' or sensor == 'L9':
            image = metadata['cloud_mask']
            band_ds = rasterio.open(image)
            band_array = band_ds.read(1)
            
            # Bit operation, source: https://www.usgs.gov/landsat-missions/landsat-collection-1-level-1-quality-assessment-band
            mask_cloud = np.bitwise_and(band_array,2**3) == 2**3
           
        # Padding 
        mask_cloud = np.pad(mask_cloud, ((0, pad_rows), (0, pad_columns)), mode='constant', constant_values=0)
            
        # Overall mask 
        masks = {}
        
        # Original resolution 
        masks[str(metadata['resolution']) + 'm'] = np.logical_or(mask_cirrus[str(metadata['resolution']) + 'm'],mask_cloud)
        
        # Add other resolution 
        for res in mask_res[1:]:
            
            # height and width at this resolution 
            height_reshaped = int( metadata['AEC_height'] / (res / metadata['resolution'])) 
            width_reshaped = int( metadata['AEC_width'] / (res / metadata['resolution'])) 
        
            # Reshape
            mask_cloud_reshaped = mask_cloud.reshape([height_reshaped, int(res / metadata['resolution']), 
                                                      width_reshaped, int(res / metadata['resolution'])]).mean(3).mean(1) > 0
            mask_reshaped = np.logical_or(mask_cirrus[str(res) + 'm'],mask_cloud_reshaped)
            masks[str(res) + 'm'] = mask_reshaped

        del mask_cirrus, mask_cloud, band_ds, band_array
        gc.collect()

        return masks
    
    # All non-water 
    elif mask_type == 'all':
    
        # SWIR band 
        print('Computing SWIR mask, reading band {}...'.format(metadata['SWIR_mask']))
        mask_SWIR = mask_threshold(metadata['SWIR_mask'],'mask_SWIR_threshold', mask_NAN = True)
        
        ### High TOA bands 
        masks_highTOA = {}
        print('Computing non-water mask, reading bands...')
        for highTOA_band_name in highTOA_band_names:
            print(highTOA_band_name)
            mask_highTOA = mask_threshold(highTOA_band_name,'mask_highTOA_threshold', mask_NAN = True)
            masks_highTOA[highTOA_band_name] = mask_highTOA
            
        # Overall mask 
        masks = {}
        
        # For each resolution 
        for res in mask_res:
            
            # Height and width at this resolution 
            height_reshaped = int( metadata['AEC_height'] / (res / metadata['resolution'])) 
            width_reshaped = int( metadata['AEC_width'] / (res / metadata['resolution'])) 
            
            # Make an empty array for each resolution 
            mask_highTOA_res = np.full((height_reshaped, width_reshaped), False)
            
            # For each highTOA band 
            for highTOA_band_name in highTOA_band_names: 
                mask_highTOA_res = np.logical_or(mask_highTOA_res, masks_highTOA[highTOA_band_name][str(res) + 'm'])
                
            # Overall mask: cirrus, SWIR, and highTOA
            mask_reshaped = np.logical_or(np.logical_or(mask_cirrus[str(res) + 'm'], mask_SWIR[str(res) + 'm']), mask_highTOA_res)
            masks[str(res) + 'm'] = mask_reshaped
        
        del mask_cirrus, mask_SWIR, mask_highTOA, masks_highTOA, mask_highTOA_res, mask_reshaped
        gc.collect()
        
        return masks   

    else:
        sys.exit('Warning: unrecognized mask was requested')

