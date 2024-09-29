# This file is part of TMart.
#
# Copyright 2024 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# Write aerosol and atmosphere information

def write_atm_info(file, basename, anci, AOT):
    import os
    
    # File path 
    path = '{}/tmart_atm_info_{}.txt'.format(file,basename)
    
    # If exist, delete
    if os.path.exists(path): 
        os.remove(path)
    
    # Write info 
    file_AER = open(path,"w")
    file_AER.write("Ratio of maritime aerosol in maritime/continental mixture: {}".format(anci['r_maritime']))
    file_AER.write("\nAerosol angstrom exponent: {}".format(anci['Angstrom_exp']))
    file_AER.write("\nAerosol single scattering albedo: {}".format(anci['SSA']))
    file_AER.write("\nAOT550: {}".format(AOT))
    file_AER.write("\nTotal column ozone: {} DU".format(anci['ozone']))
    file_AER.write("\nTotal precipitable water vapour: {} kg m-2".format(anci['water_vapour']))
    file_AER.close()
