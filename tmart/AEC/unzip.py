# This file is part of TMart.
#
# Copyright 2024 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Unzip a .SAFE.zip file and return the path to the extracted .SAFE dir
# With contribution from Shun Bi

def unzip(zip_file_path):
    import zipfile
    unzip_folder = zip_file_path.replace('.SAFE.zip', '.SAFE')
    with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
        zip_ref.extractall(unzip_folder)
    print(f"Unzipped {zip_file_path} to {unzip_folder}")
    return unzip_folder
