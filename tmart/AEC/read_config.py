# This file is part of TMart.
#
# Copyright 2024 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Read configuration file 

def read_config(mask_SWIR_threshold):
    import os
    import pandas as pd
    
    # Directory 
    one_up = os.path.dirname(os.path.dirname(__file__ ))
    file = os.path.join(one_up, 'config/config.txt')

    # Read file 
    df = pd.read_csv(file, skiprows=0, skipinitialspace=True,on_bad_lines='skip')

    # Remvoe comments 
    df = df.drop(df[df.iloc[:,0].str.startswith('#')].index)
    
    # Extract info 
    df[['left', 'right']] = df.iloc[:,0].str.replace(' ', '').str.split('=', expand=True)
    result_dict = dict(zip(df['left'], df['right']))
    print('\nT-Mart configuration (editable at {}): '.format(file))
    
    if mask_SWIR_threshold is not None:
        result_dict['mask_SWIR_threshold']=mask_SWIR_threshold
    
    for k, v in result_dict.items():
        print(str(k) + ': '  + str(v))    
    
    return result_dict
    