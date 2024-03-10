# This file is part of T-Mart.
#
# Copyright 2023 Yulun Wu.
#
# T-Mart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import numpy as np

# Calculate the absorption optical thickness between two points 
def find_OT(q0,q1,atm_profile):
 
    topP = max(q0[2],q1[2])
    bottomP = min(q0[2],q1[2])
    
    # Full layers to include 
    top_below_topP = atm_profile[:,1] <= topP/1000
    bottom_above_bottomP = atm_profile[:,0] >= bottomP/1000 
    full_layers = top_below_topP & bottom_above_bottomP
    tao_abs = sum(atm_profile[full_layers,2])

    # TopP 
    bottom_below_topP = atm_profile[:,0] < topP/1000 # bottom of layer is lower than z 
    top_above_topP = atm_profile[:,1] > topP/1000 # top of layer is above z 
    within_topP = bottom_below_topP & top_above_topP
    
    # BottomP 
    bottom_below_bottomP = atm_profile[:,0] < bottomP/1000 
    top_above_bottomP = atm_profile[:,1] > bottomP/1000 
    within_bottomP = bottom_below_bottomP & top_above_bottomP
    
    # Same layer 
    if np.any(within_topP) and np.any(within_bottomP) and np.all(within_topP == within_bottomP):
        
        remain_ratio = (topP-bottomP)/1000 / atm_profile[within_topP,6][0]
        tao_abs = atm_profile[within_topP,2][0] * remain_ratio
    
    else:
        # Top remain
        if np.any(within_topP):
            remain_ratio = (topP/1000 - atm_profile[within_topP,0][0]) / atm_profile[within_topP,6][0] # remain / total height
            tao_abs = tao_abs + atm_profile[within_topP,2][0] * remain_ratio
            
        # Bottom remain 
        if np.any(within_bottomP):
            remain_ratio = (atm_profile[within_bottomP,1][0] - bottomP/1000  ) / atm_profile[within_bottomP,6][0]
            tao_abs = tao_abs + atm_profile[within_bottomP,2][0] * remain_ratio
        
    return tao_abs
