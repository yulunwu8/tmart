# This file is part of TMart.
#
# Copyright 2024-2025 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# Read aerosol and atmosphere information from a text file

def read_atm_info(path):
    import os
    import re

    if not os.path.exists(path):
        raise FileNotFoundError(f'Atmospheric info file not found: {path}')

    label_map = {
        'Ratio of maritime aerosol in maritime/continental mixture': 'r_maritime',
        'Aerosol angstrom exponent': 'Angstrom_exp',
        'Aerosol single scattering albedo': 'SSA',
        'AOT550': 'AOT_MERRA2',
        'Total column ozone': 'ozone',
        'Total precipitable water vapour': 'water_vapour',
    }

    values = {key: None for key in label_map.values()}

    with open(path, 'r') as atm_file:
        for line in atm_file:
            for label, key in label_map.items():
                if label in line:
                    body = line.split(':', 1)[1]
                    match = re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', body)
                    if not match:
                        raise ValueError(f'No numeric value found for "{label}" in {path}')
                    values[key] = float(match[0])

    required_keys = ['r_maritime', 'ozone', 'water_vapour']
    missing = [key for key in required_keys if values.get(key) is None]
    if missing:
        missing_list = ', '.join(missing)
        raise ValueError(f'Missing required fields in {path}: {missing_list}')

    return values
