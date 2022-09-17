# This file is part of TMart.
#
# Copyright 2022 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.





# Aerosol SPF

import numpy as np
import pandas as pd
import os.path

def find_aerosolSPF(aerosol_type,wl):

    # wl = 3750
    
    # Currently only supporting the default number of angles in 6S
    # For higher angular resolution SPF or mixed aerosols, 
    # contact Yulun Wu at yulunwu8@gmail.com
    
    # Currently SPF is interpolated using the central wavelength. Alternatively,
    # can run 6S and extract the spectually resampled values. 
    
    
    n_angles = 83
    
    ### Columns 
    
    # From 6S code and manual 
    wls = np.array([0.350,0.400,0.412,0.443,0.470,0.488,0.515,0.550,0.590,0.633,0.670,0.694,0.760,0.860,1.240,1.536,1.650,1.950,2.250,3.750]) * 1000
    
    ### Rows, the angles
    
    angles = np.empty(83, dtype=float)
    
    angles[0] = 0.
    angles[1:41] = np.arange(0.75, 40, 1) / 80.5 * 180
    angles[41]= 90.
    angles[42:82] = np.arange(40.25 + 0.5, 80, 1) / 80.5 * 180 # 40.25 is right in the middle
    angles[82]=180.
    
    # reverse the angles 
    angles[::-1].sort()
    
    ###
    

    file_aerosolSPF = 'ancillary/aerosolSPF/' + str(aerosol_type) + '.csv'
    file_aerosolSPF = os.path.join(os.path.dirname(__file__), file_aerosolSPF)
    
    aerosolSPF = np.genfromtxt(file_aerosolSPF, delimiter=',')
    
    # np.shape(aerosolSPF)
    # aerosolSPF[0,:]
    # np.interp(wl, wls, aerosolSPF[0,:])
    
    
    aerosolSPF_wl = np.array([np.interp(wl, wls, aerosolSPF[i,:]) for i in range(np.shape(aerosolSPF)[0])])
    aerosolSPF_wl = aerosolSPF_wl[0:n_angles]
    
    
    df = pd.DataFrame({'Angle':angles, 'Value':aerosolSPF_wl}).sort_values('Angle').reset_index()
    
    # df
    
    # df.plot( 'angle','value', kind='scatter', logy=True)
    
    
    # Already normalized data
    
    
    # df.to_csv('testSPF.csv', index=False)
    return df


if __name__=='__main__':
    aerosol_type = 'Maritime'
    # aerosol_type = 'Continental'
    
    wl = 550
    
    test = find_aerosolSPF(aerosol_type,wl)
    print(test)







