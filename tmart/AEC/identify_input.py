# This file is part of TMart.
#
# Copyright 2024 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Identify the type of the input: directory, file, or .SAFE.zip
# With contribution from Shun Bi

def identify_input(file_path):
    import os, sys
    import tmart
    
    # Initialize file_is_dir
    file_is_dir = None
    
    # If input is a directory
    if os.path.isdir(file_path):
        file_is_dir = True
        
    # If input is a file 
    elif os.path.isfile(file_path):
        # ACOLITE L1R (currently supporting PRISMA data only)
        if file_path.endswith('L1R.nc'):
            file_is_dir = False
        # SAFE.zip file: unzip and update the path 
        elif file_path.endswith('.SAFE.zip'):
            file_path = tmart.AEC.unzip(file_path)
            file_is_dir = True
        
    # If the file or directory was not identified
    if file_is_dir is None:
        sys.exit(f"File or folder cannot be identified: {file_path}")
        
    return file_path, file_is_dir
