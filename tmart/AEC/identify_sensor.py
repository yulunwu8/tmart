# This file is part of TMart.
#
# Copyright 2024 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# Identify satellite sensor, currently only supports S2 and L8/L9 

def identify_sensor(file):
    import os, sys
    sensor = None

    try: 
        strings = file.split(".")
        base_name = os.path.basename(file)
        
        # remove AEC_ if exists 
        if base_name[0:4] == 'AEC_':
            base_name = base_name[4:]

        # if with 'SAFE': S2
        if strings[-1] == 'SAFE':
            if   base_name[0:3] == 'S2A': sensor = 'S2A'
            elif base_name[0:3] == 'S2B': sensor = 'S2B'
            elif base_name[0:3] == 'S2C': sensor = 'S2C'
        
        # else: Landsat series 
        else:
            if   base_name[0:4] == 'LC08': sensor = 'L8'
            elif base_name[0:4] == 'LC09': sensor = 'L9'

    except:
        pass

    if sensor is None:
        sys.exit('Unable to identify sensor, please use the original folder name.')
    
    print('\nIdentified sensor: ' + str(sensor) +  '\n')

    return sensor
    
    