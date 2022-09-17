# This file is part of TMart.
#
# Copyright 2022 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.







# Photon movement


import numpy as np
import math
import sys

from .tm_geometry import dirC_to_coord, dirP_to_coord
      
        
def pt_move(atm_profile,q0,pt_direction,sampled_tao):        
    '''
    

    Parameters
    ----------
    atm_profile : TYPE
        Scattering and absorptions coefficients in the atmosphere.
    q0 : TYPE
        Starting point.
    pt_direction : TYPE
        The moving direction.
    sampled_tao : TYPE
        The optical thickness to move in space.

    Returns
    -------
    q1 : TYPE
        end point of the movement .
    out : Boolean
        If out of atmosphere after the movement.

    '''
        
    # If the photon is out of the atmosphere after the movement 
    out = False
    n_layers = len(atm_profile)
    
    # Percentage traveled in each layer 
    atm_profile[:,7] = 0
    
    if pt_direction[0]<0 or pt_direction[0] > 180:
        print('pt_direction: ' + str(pt_direction))
        sys.exit('pt_direction error')  
    
    ### In case zenith = 90, this should be impossible but we'll process it 
    elif pt_direction[0]==90:
        
        print ('Photon travel exactly parallel to the atm layers, check absorption...')

        # Identify layers 
        bottom_below = atm_profile[:,0] < q0[2]/1000 # bottom of layer is lower than z 
        top_above = atm_profile[:,1] > q0[2]/1000 # top of layer is above z 
        within = bottom_below & top_above       
        
        # layer the photon is in 
        n_within = np.where(within)[0][0] 
        c = atm_profile[n_within,5] / atm_profile[n_within,6]  # ot_scatt / layer_height in km
        
        traveled_distance = sampled_tao / c    * 1000  # convert to m 
        
        q1 = dirC_to_coord (dirP_to_coord(1, pt_direction),q0,traveled_distance)
        
        tao_abs =       sampled_tao * ( atm_profile[n_within,2] / atm_profile[n_within,5] ) # ot_abs / layer_height 
        ot_rayleigh =   sampled_tao * ( atm_profile[n_within,3] / atm_profile[n_within,5] ) 
        ot_mie =        sampled_tao * ( atm_profile[n_within,4] / atm_profile[n_within,5] ) 
        
        return q1, tao_abs, ot_rayleigh, ot_mie, out 
        
    
    ### Move down if > 90
    elif pt_direction[0]>90:
    
        travel_angle = 180 - pt_direction[0] # relative to parallel plane
        
        # Identify layers 
        top_below = atm_profile[:,1] <= q0[2]/1000  # top of layer is equal or lower than z
        bottom_below = atm_profile[:,0] < q0[2]/1000 # bottom of layer is lower than z 
        top_above = atm_profile[:,1] > q0[2]/1000 # top of layer is above z 
        within = bottom_below & top_above  
    
        # If movement within a layer
        if np.any(within):
            
            # Layer the photon is in 
            n_within = np.where(within)[0][0] 
            
            # The bottom of the layer within 
            within_bottom = atm_profile[within,0][0] 
 
            # leftover_height / layer height
            leftover_height_ratio = (q0[2]/1000 - within_bottom) /  atm_profile[within,6][0] 
       
            atm_profile[within,7] = leftover_height_ratio / math.cos(travel_angle/180*math.pi)
        
        atm_profile[top_below,7] = 1 / math.cos(travel_angle/180*math.pi)      
        
        # Sum of ot_scatt * percentage 
        tao_sum = sum(atm_profile[:,5]*atm_profile[:,7]) 

        # If traveled tao is larger than what's there, AKA crossing layers 
        if sampled_tao > tao_sum:
            z = -10 # penetrate to 10m underground 
            tao_abs =       sum(atm_profile[:,2]*atm_profile[:,7])
            ot_rayleigh =   sum(atm_profile[:,3]*atm_profile[:,7])
            ot_mie =        sum(atm_profile[:,4]*atm_profile[:,7])
            
        # If movement happened within the layer 
        else: 
            
            ## Remove layers one by one until we reach the sampled_tao
            for i in range(0,n_layers): # remove layer i
            
                # From i+1 to last layer  
                atm_profile_temp = atm_profile[ range(i+1,n_layers)  ,:]   
                
                # Recalculate tao_sum
                tao_sum = sum(atm_profile_temp[:,5]*atm_profile_temp[:,7])
                
                # Once we find the layer 
                if sampled_tao > tao_sum:
                    
                    tao_lastLayer = sampled_tao - tao_sum  
                    
                    # If movement within a layer
                    if np.any(within):
                        if i  == n_within: 
                            tao_lastLayer_ratio = tao_lastLayer / (atm_profile[i,5] / math.cos(travel_angle/180*math.pi) ) 
                            
                            z = q0[2] - atm_profile[i,6] * tao_lastLayer_ratio * 1000
                            tao_abs =       atm_profile[i,2] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)
                            ot_rayleigh =   atm_profile[i,3] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)
                            ot_mie =        atm_profile[i,4] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)
                            
                            break
                    
                    tao_lastLayer_ratio = tao_lastLayer / (atm_profile[i,5] * atm_profile[i,7]) 
                    z = atm_profile[i,1] - atm_profile[i,6] * tao_lastLayer_ratio # top of the layer - proportion in the last layer 
                    z = z * 1000
                    
                    tao_abs =       sum(atm_profile_temp[:,2]*atm_profile_temp[:,7]) + atm_profile[i,2] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)
                    ot_rayleigh =   sum(atm_profile_temp[:,3]*atm_profile_temp[:,7]) + atm_profile[i,3] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)
                    ot_mie =        sum(atm_profile_temp[:,4]*atm_profile_temp[:,7]) + atm_profile[i,4] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)
                    
                    break # stop removing layers 
    
    ### Move up if < 90
    elif pt_direction[0]<90: 
    
        travel_angle = pt_direction[0] # relative to parallel plane
        
        # Identify layers 
        bottom_above = atm_profile[:,0] >= q0[2]/1000 # bottom of layer is equal or above z     
        bottom_below = atm_profile[:,0] < q0[2]/1000 # bottom of layer is lower than z 
        top_above = atm_profile[:,1] > q0[2]/1000 # top of layer is above z 
        within = bottom_below & top_above  

        # If movement within a layer
        if np.any(within):

            # Layer the photon is in 
            n_within = np.where(within)[0][0]         
            
            # The top of the layer within 
            within_top = atm_profile[within,1][0] 
            
            # leftover_height / layer height
            leftover_height_ratio = (within_top - q0[2]/1000 ) /  atm_profile[within,6][0] 
            
            atm_profile[within,7] = leftover_height_ratio / math.cos(travel_angle/180*math.pi)
            
        atm_profile[bottom_above,7] = 1 / math.cos(travel_angle/180*math.pi)
      
        # Sum of ot_scatt * percentage 
        tao_sum = sum(atm_profile[:,5]*atm_profile[:,7]) 
        
        # If traveled tao is larger than what's there, AKA crossing layers 
        if sampled_tao > tao_sum:
            z = max(atm_profile[:,1]) * 1000  # penetrate through TOA 
            tao_abs =       sum(atm_profile[:,2]*atm_profile[:,7])
            ot_rayleigh =   sum(atm_profile[:,3]*atm_profile[:,7])
            ot_mie =        sum(atm_profile[:,4]*atm_profile[:,7]) 
            out = True
             
        # If movement happened within the layer    
        else: 
        
            ## Remove layers one by one until we reach the sampled_tao
            for i in range(0,n_layers): 

                # From i+1 to last layer   
                atm_profile_temp = atm_profile[ range(0,n_layers-1-i) ,:]  
                
                # Recalculate tao_sum
                tao_sum = sum(atm_profile_temp[:,5]*atm_profile_temp[:,7])
                
                # Once we find the layer
                if sampled_tao > tao_sum:

                    layer_index = n_layers-i-1
                    tao_lastLayer = sampled_tao - tao_sum 
                    
                    # If movement within a layer
                    if np.any(within):
                        if layer_index  == n_within: 
                            tao_lastLayer_ratio = tao_lastLayer / (atm_profile[layer_index,5] / math.cos(travel_angle/180*math.pi) ) 
                            
                            z = q0[2] + atm_profile[layer_index,6] * tao_lastLayer_ratio * 1000
                            tao_abs =       atm_profile[layer_index,2] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)
                            ot_rayleigh =   atm_profile[layer_index,3] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)
                            ot_mie =        atm_profile[layer_index,4] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)
                            
                            break
                        
                    tao_lastLayer_ratio = tao_lastLayer / (atm_profile[layer_index,5] * atm_profile[layer_index,7]) 
                    z = atm_profile[layer_index,0] + atm_profile[layer_index,6] * tao_lastLayer_ratio # bottom of the layer + proportion in the last layer 
                    z = z * 1000
                    
                    tao_abs =       sum(atm_profile_temp[:,2]*atm_profile_temp[:,7]) + atm_profile[layer_index,2] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)   
                    ot_rayleigh =   sum(atm_profile_temp[:,3]*atm_profile_temp[:,7]) + atm_profile[layer_index,3] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)
                    ot_mie =        sum(atm_profile_temp[:,4]*atm_profile_temp[:,7]) + atm_profile[layer_index,4] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)
                    
                    break # stop removing layers 
           
    traveled_distance = abs(q0[2] - z) / math.cos(travel_angle/180*math.pi)
    q1 = dirC_to_coord (dirP_to_coord(1, pt_direction),q0,traveled_distance)
    q1 = np.array(q1)
    
    return q1, tao_abs, ot_rayleigh, ot_mie, out
        





# if __name__=='__main__':  
    
#     atm_profile = my_tmart.atm_profile_wl.sort_values('Alt_bottom').to_numpy()
    
#     q0 = np.array([0.0, 0.0, 90_000.0])
#     pt_direction = [180, 0]
#     sampled_tao = 1e-5
    
#     q1, tao_abs, ot_rayleigh, ot_mie, out = pt_move(atm_profile,q0,pt_direction,sampled_tao)
    
#     print('q1: ' + str(q1)) 
#     print('tao_abs: ' + str(tao_abs)) 
#     print('ot_rayleigh: ' + str(ot_rayleigh)) 
#     print('ot_mie: ' + str(ot_mie)) 
#     print('out: ' + str(out)) 
    
    
    
    
    
    
        


