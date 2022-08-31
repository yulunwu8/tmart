


import pandas as pd
import numpy as np


# calculate the optical thickness between two points 

def find_OT(q0,q1,atm_profile):
 
    
    topP = max(q0[2],q1[2])
    bottomP = min(q0[2],q1[2])
    
    # full layers to include 
    
    top_below_topP = atm_profile[:,1] <= topP/1000
    # print('top_below_topP: ' + str(top_below_topP))
    
    bottom_above_bottomP = atm_profile[:,0] >= bottomP/1000 
    # print('bottom_above_bottomP: ' + str(bottom_above_bottomP))
    
    
    full_layers = top_below_topP & bottom_above_bottomP
    # print('full_layers: ' + str(full_layers))
    
    tao_abs = sum(atm_profile[full_layers,2])
    # print('tao_abs: ' + str(tao_abs))
    
    
    # topP 
    
    bottom_below_topP = atm_profile[:,0] < topP/1000 # bottom of layer is lower than z 
    top_above_topP = atm_profile[:,1] > topP/1000 # top of layer is above z 
    within_topP = bottom_below_topP & top_above_topP
    
    # print('\ntopP: ' + str(topP))
    # print('within_topP: ' + str(within_topP))
    
    
    # bottomP 
    
    bottom_below_bottomP = atm_profile[:,0] < bottomP/1000 
    top_above_bottomP = atm_profile[:,1] > bottomP/1000 
    within_bottomP = bottom_below_bottomP & top_above_bottomP
    
    # print('bottomP: ' + str(bottomP))
    # print('within_bottomP: ' + str(within_bottomP))
    
    
    
    if np.any(within_topP) and np.any(within_bottomP) and np.all(within_topP == within_bottomP):
        # print('\nSame layer')
        
        remain_ratio = (topP-bottomP)/1000 / atm_profile[within_topP,6][0]
        # print('remain_ratio: ' + str(remain_ratio))
        
        tao_abs = atm_profile[within_topP,2][0] * remain_ratio
        # print('tao_abs: ' + str(tao_abs))
          
    
    else:
        
        if np.any(within_topP):
            # print('\nTop remain')
            
            remain_ratio = (topP/1000 - atm_profile[within_topP,0][0]) / atm_profile[within_topP,6][0] # remain / total height
            # print('remain_ratio: ' + str(remain_ratio))
            
            tao_abs = tao_abs + atm_profile[within_topP,2][0] * remain_ratio
            # print('tao_abs: ' + str(tao_abs))
            
        if np.any(within_bottomP):
            # print('\nBottom remain')
            
            remain_ratio = (atm_profile[within_bottomP,1][0] - bottomP/1000  ) / atm_profile[within_bottomP,6][0]
            # print('remain_ratio: ' + str(remain_ratio))
            
            tao_abs = tao_abs + atm_profile[within_bottomP,2][0] * remain_ratio
            # print('tao_abs: ' + str(tao_abs))
        
    
    return tao_abs




if __name__=='__main__':

    
    atm_profile_wl = pd.read_csv('test_atm_profile.csv')
    atm_profile_wl = pd.read_csv('Atm_Rayleigh014_test.csv')
    
    
    
    atm_profile_wl['ot_scatt'] = atm_profile_wl.ot_rayleigh + atm_profile_wl.ot_mie
    atm_profile_wl['l_height'] = atm_profile_wl.Alt_top - atm_profile_wl.Alt_bottom
    atm_profile_wl['percentage'] = 0 # used to calculate travelling distance 
    
    atm_profile = atm_profile_wl.sort_values('Alt_bottom').to_numpy()
    

    q0 = [0,0,0]
    
    q1 = [100,100, 3_000]
    
    tao_abs = find_OT(q0,q1,atm_profile)
    print(tao_abs)
    
    


























