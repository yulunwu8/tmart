# This file is part of TMart.
#
# Copyright 2025 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# Plot the extent of water in the satellite image folder

def plot_water_extent(metadata, config, mask_all):
    
    import rasterio, sys, os, math, re, gc
    import numpy as np
    import matplotlib.pyplot as plt
    
    # Base name of sat. file
    base_name = os.path.basename(metadata['file'])
    base_name = re.sub(r'^AEC_', '', base_name)
    base_name = re.sub(r'\.SAFE$', '', base_name)
    
    # Identify bands
    sensor = metadata['sensor']
    if sensor == 'S2A' or sensor == 'S2B' or sensor == 'S2C':
        RGB = ['B04','B03','B02']
        res = '10m'
    elif sensor == 'L8' or sensor == 'L9':
        RGB = ['B4','B3','B2']
        res = '30m'
    
    image_bands = []    
    
    # Loop through 3 colours 
    for band_name in RGB:        
        band_file = metadata[band_name]
        band_ds = rasterio.open(band_file)
        
        # Scaling  
        scale_mult = metadata[str(band_name) + '_mult']
        scale_add = metadata[str(band_name) + '_add']
        
        
        # Read image
        band_ds_image = band_ds.read(1)
        image = band_ds_image * scale_mult + scale_add
        
        # L8 solar zenith correction: https://www.usgs.gov/landsat-missions/using-usgs-landsat-level-1-data-product
        if sensor =='L8' or sensor == 'L9': image = image / math.cos(metadata['sza']/180*math.pi)
    
        image_bands.append(image)
    
    # Load mask 
    mask = mask_all[res]
    if mask.shape != image.shape:
        mask = mask[:image.shape[0], :image.shape[1]]
    
    # Identify % water
    is_nan = band_ds_image==0
    percent_water = np.sum(mask == False) / np.sum(is_nan == False) * 100
    print(f"{percent_water:.2f}% of non-NA pixels are identified as water.")
    
    # Stack rgb
    image_rgb = np.stack(image_bands, axis=-1)  # Stack into (H, W, 3) format
    image_rgb = np.clip(image_rgb / 0.25, 0, 1)  # Normalize and clip values
    
    # Create an overlay image 
    overlay = np.zeros_like(image_rgb)
    overlay[..., 1] = 1  # full intensity green
    overlay[..., 2] = 1  # full intensity blue
    
    # Create a per-pixel alpha channel
    alpha_channel = np.where(mask, 0.0, 1)
    
    # Combine overlay and alpha_channel to form an RGBA image (H, W, 4)
    overlaya = np.dstack((overlay, alpha_channel))
    
    # Plot
    fig, ax = plt.subplots(figsize=(8, 9))
    fig.subplots_adjust(left=0.05, right=0.95, top=0.97, bottom=0.07)
    ax.imshow(image_rgb)
    ax.imshow(overlaya)  
    ax.axis('off')       
    
    # Add annotation 
    fig.text(0.5, 0.95, base_name, weight='bold', ha='center', fontsize=12, wrap=True)
    fig.text(0.5, 0.045, 'Only pixels shaded in cyan were identified as water and were corrected for the adjacency effect by T-Mart. See <T-Mart user guide>/<Adjacency-Effect Correction>/<AEC configuration> for more information. Current mask_SWIR_threshold: ' + str(config['mask_SWIR_threshold']),
             ha='center', fontsize=12, wrap=True)
    
    # Save the figure 
    path_preview = os.path.join(metadata['file'],'tmart_preview.png')
    plt.savefig(path_preview, dpi=300, pad_inches=0.1)
    plt.close(fig)
    
    print('Preview saved: ' + str(path_preview))

    del image_bands, mask, image_rgb, overlay, alpha_channel, overlaya
    gc.collect()
