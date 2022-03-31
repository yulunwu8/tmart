# Photon movement



import numpy as np
import random, math
import sys

from .tm_geometry import dirC_to_coord, dirP_to_coord
      
        
def pt_move(atm_profile,q0,pt_direction,sampled_tao):        
    '''
    

    Parameters
    ----------
    atm_profile : TYPE
        DESCRIPTION.
    q0 : TYPE
        DESCRIPTION.
    pt_direction : TYPE
        DESCRIPTION.
    sampled_tao : TYPE
        DESCRIPTION.

    Returns
    -------
    q1 : TYPE
        end point of the movement .
    out : Boolean
        If out of atm..

    '''
        
    out = False
    
    
    ### Every calculating distance 
    n_layers = len(atm_profile)
    atm_profile[:,7] = 0
    
    # print('\natm_profile: ' + str(atm_profile))
    

    # top_below = np.array( self.atm_profile_wl.Alt_top <= q0[2]/1000 )  
    
    # print('top_below: ' + str(top_below))
    
    # print('sampled_tao: ' + str(sampled_tao))
    
    if pt_direction[0]<0 or pt_direction[0] > 180:
        print('pt_direction: ' + str(pt_direction))
        sys.exit('pt_direction error')  
    
    
    ### In case zenith = 90, this should be impossible 
    elif pt_direction[0]==90:
        
        print ('Photon travel exactly parallel to the atm layers, check absorption...')

        
        bottom_below = atm_profile[:,0] < q0[2]/1000 # bottom of layer is lower than z 
        top_above = atm_profile[:,1] > q0[2]/1000 # top of layer is above z 
        within = bottom_below & top_above       
        
        # print('within: ' + str(within))
        
        n_within = np.where(within)[0][0] # layer the photon is in 
        c = atm_profile[n_within,5] / atm_profile[n_within,6]  # ot_scatt / layer_height 
        
        
        traveled_distance = sampled_tao / c    * 1000 
        
        # print('traveled_distance: ' + str(traveled_distance/20000))
        
        q1 = dirC_to_coord (dirP_to_coord(1, pt_direction),q0,traveled_distance)
        
        tao_abs =       sampled_tao * ( atm_profile[n_within,2] / atm_profile[n_within,5] ) # ot_abs / layer_height 
        ot_rayleigh =   sampled_tao * ( atm_profile[n_within,3] / atm_profile[n_within,5] ) 
        ot_mie =        sampled_tao * ( atm_profile[n_within,4] / atm_profile[n_within,5] ) 
        
        
        return q1, tao_abs, ot_rayleigh, ot_mie, out 
        
        
    
    
    ### move down if > 90
    elif pt_direction[0]>90:
    
        travel_angle = 180 - pt_direction[0] # relative to parallel plane
        # print('travel_angle: ' + str(travel_angle))
        
        
        
        top_below = atm_profile[:,1] <= q0[2]/1000  # top of layer is equal or lower than z
        
        bottom_below = atm_profile[:,0] < q0[2]/1000 # bottom of layer is lower than z 
        top_above = atm_profile[:,1] > q0[2]/1000 # top of layer is above z 
        within = bottom_below & top_above  
        
        # print('\n')
        

        # print('top_below: ' + str(top_below))
        # print('bottom_below: ' + str(bottom_below))
        # print('top_above: ' + str(top_above))
        # print('within: ' + str(within))
        

        # if within a layer, add this 
        if np.any(within):
            # print('\nyes within')
            
            n_within = np.where(within)[0][0] # layer the photon is in 
            # print('n_within: ' + str(n_within))
            
            within_bottom = atm_profile[within,0][0] # the bottom of the layer within 
            # print('within_bottom: ' + str(within_bottom))
            
            leftover_height_ratio = (q0[2]/1000 - within_bottom) /  atm_profile[within,6][0] # leftover_height / layer height
            # print('leftover_height_ratio: ' + str(leftover_height_ratio))
            
            atm_profile[within,7] = leftover_height_ratio / math.cos(travel_angle/180*math.pi)
        

        atm_profile[top_below,7] = 1 / math.cos(travel_angle/180*math.pi)
        # print('atm_profile: ' + str(atm_profile))        
        
        
        tao_sum = sum(atm_profile[:,5]*atm_profile[:,7]) # sum of ot_scatt * percentage 
        # print('tao_sum: ' + str(tao_sum))
        
        
        if sampled_tao > tao_sum:
            
            z = -10 # penetrate to 10m underground 
            
            tao_abs =       sum(atm_profile[:,2]*atm_profile[:,7])
            # print('tao_abs: ' + str(tao_abs))
            ot_rayleigh =   sum(atm_profile[:,3]*atm_profile[:,7])
            ot_mie =        sum(atm_profile[:,4]*atm_profile[:,7])
            
            
            
        else: # within layer 
        
            # print ('\nsampled_tao < tao_sum')
            
            ## for loop, number of layers in the atm 
            
            for i in range(0,n_layers): # remove layer i
                # print('====== ' + str(i))
                
                atm_profile_temp = atm_profile[ range(i+1,n_layers)  ,:] # from i+1 to last layer    
                # print('atm_profile_temp: ' + str(atm_profile_temp))
                
                # recalculate tao_sum
                tao_sum = sum(atm_profile_temp[:,5]*atm_profile_temp[:,7])
                # print('tao_sum: ' + str(tao_sum)) 
                
                
                if sampled_tao > tao_sum:
                    # print('\nThis is the layer')
                    
                    tao_lastLayer = sampled_tao - tao_sum 
                    # print('tao_lastLayer: ' + str(tao_lastLayer))     

                
                    if np.any(within):
                        if i  == n_within: # only travel within layer
                            # print('\ntravel in the same layer')
                            tao_lastLayer_ratio = tao_lastLayer / (atm_profile[i,5] / math.cos(travel_angle/180*math.pi) ) 
                            # print('tao_lastLayer_ratio: ' + str(tao_lastLayer_ratio)) 
                            
                            z = q0[2] - atm_profile[i,6] * tao_lastLayer_ratio * 1000
                            
                            tao_abs =       atm_profile[i,2] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)
                            # print('tao_abs: ' + str(tao_abs)) 
                            
                            ot_rayleigh =   atm_profile[i,3] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)
                            ot_mie =        atm_profile[i,4] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)
                            
                            break
                    

                    tao_lastLayer_ratio = tao_lastLayer / (atm_profile[i,5] * atm_profile[i,7]) # proportion in the last layer 
                    # print('tao_lastLayer_ratio: ' + str(tao_lastLayer_ratio)) 
                    
                    
                    z = atm_profile[i,1] - atm_profile[i,6] * tao_lastLayer_ratio # top of the layer - proportion in the last layer 
                    z = z * 1000
                    
                    
                    tao_abs =       sum(atm_profile_temp[:,2]*atm_profile_temp[:,7]) + atm_profile[i,2] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)
                    # print('tao_abs: ' + str(tao_abs)) 
                    ot_rayleigh =   sum(atm_profile_temp[:,3]*atm_profile_temp[:,7]) + atm_profile[i,3] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)
                    ot_mie =        sum(atm_profile_temp[:,4]*atm_profile_temp[:,7]) + atm_profile[i,4] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)
                    
                    break # stop removing layers 
    
    # move up 
    elif pt_direction[0]<90: 
    
        travel_angle = pt_direction[0] # relative to parallel plane
        # print('travel_angle: ' + str(travel_angle))

        
        bottom_above = atm_profile[:,0] >= q0[2]/1000 # bottom of layer is equal or above z
        
        bottom_below = atm_profile[:,0] < q0[2]/1000 # bottom of layer is lower than z 
        top_above = atm_profile[:,1] > q0[2]/1000 # top of layer is above z 
        within = bottom_below & top_above  
        
        # print('\n')
        # print('bottom_above: ' + str(bottom_above))
        # print('bottom_below: ' + str(bottom_below))
        # print('top_above: ' + str(top_above))
        # print('within: ' + str(within))
        

        # if within a layer, add this 
        if np.any(within):
            # print('\nyes within')
            
            n_within = np.where(within)[0][0] # layer the photon is in 
            # print('n_within: ' + str(n_within))            
            
            within_top = atm_profile[within,1][0] # the top of the layer the photon is in  
            # print('within_top: ' + str(within_top))
            
            leftover_height_ratio = (within_top - q0[2]/1000 ) /  atm_profile[within,6][0] # leftover_height / layer height
            # print('leftover_height_ratio: ' + str(leftover_height_ratio))
            
            atm_profile[within,7] = leftover_height_ratio / math.cos(travel_angle/180*math.pi)
            

        
        atm_profile[bottom_above,7] = 1 / math.cos(travel_angle/180*math.pi)
        # print('atm_profile: ' + str(atm_profile))        
        
        
        tao_sum = sum(atm_profile[:,5]*atm_profile[:,7]) # sum of ot_scatt * percentage 
        # print('tao_sum: ' + str(tao_sum))
        
        
        # if traveled tao is larger than what's there 
        if sampled_tao > tao_sum:
            z = max(atm_profile[:,1]) * 1000  # penetrate through TOA 
            # print('z: ' + str(z)) 
            
            out = True
            
            tao_abs =       sum(atm_profile[:,2]*atm_profile[:,7])
            # print('tao_abs: ' + str(tao_abs))
            ot_rayleigh =   sum(atm_profile[:,3]*atm_profile[:,7])
            ot_mie =        sum(atm_profile[:,4]*atm_profile[:,7])            
            

         
            
        else: # within layer 
        
            # print ('\nsampled_tao < tao_sum')
            
            ## for loop, number of layers in the atm 
            
            for i in range(0,n_layers): # remove layer i
                # print('====== ' + str(i))
                
                atm_profile_temp = atm_profile[ range(0,n_layers-1-i) ,:] # from i+1 to last layer    
                # print('atm_profile_temp: ' + str(atm_profile_temp))
                
                # recalculate tao_sum
                tao_sum = sum(atm_profile_temp[:,5]*atm_profile_temp[:,7])
                # print('tao_sum: ' + str(tao_sum)) 
    
                if sampled_tao > tao_sum:
                    # print('\nThis is the layer')
                    
                    layer_index = n_layers-i-1
                    # print('layer_index: ' + str(layer_index)) 
                    
                    tao_lastLayer = sampled_tao - tao_sum 
                    # print('tao_lastLayer: ' + str(tao_lastLayer)) 
                    
                    # print('======')
                    # print(atm_profile[n_layers-i,:])
                    
                    
                    if np.any(within):
                        if layer_index  == n_within: # only travel within layer
                            # print('\ntravel in the same layer')
                            tao_lastLayer_ratio = tao_lastLayer / (atm_profile[layer_index,5] / math.cos(travel_angle/180*math.pi) ) 
                            # print('tao_lastLayer_ratio: ' + str(tao_lastLayer_ratio)) 
                            
                            z = q0[2] + atm_profile[layer_index,6] * tao_lastLayer_ratio * 1000
                            
                            tao_abs =       atm_profile[layer_index,2] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)
                            # print('tao_abs: ' + str(tao_abs)) 
                            ot_rayleigh =   atm_profile[layer_index,3] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)
                            ot_mie =        atm_profile[layer_index,4] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)
                            
                            break
                        
                    
                    tao_lastLayer_ratio = tao_lastLayer / (atm_profile[layer_index,5] * atm_profile[layer_index,7]) # proportion in the last layer 
                    # print('tao_lastLayer_ratio: ' + str(tao_lastLayer_ratio)) 
                    
                    
                    z = atm_profile[layer_index,0] + atm_profile[layer_index,6] * tao_lastLayer_ratio # bottom of the layer + proportion in the last layer 
                    z = z * 1000
                    
                    tao_abs =       sum(atm_profile_temp[:,2]*atm_profile_temp[:,7]) + atm_profile[layer_index,2] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)
                    # print('tao_abs: ' + str(tao_abs))      
                    ot_rayleigh =   sum(atm_profile_temp[:,3]*atm_profile_temp[:,7]) + atm_profile[layer_index,3] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)
                    ot_mie =        sum(atm_profile_temp[:,4]*atm_profile_temp[:,7]) + atm_profile[layer_index,4] * tao_lastLayer_ratio / math.cos(travel_angle/180*math.pi)
                    
                    break # stop removing layers 
    
                
    # print('\n')
            
    # print('z: ' + str(z)) 
           
    traveled_distance = abs(q0[2] - z) / math.cos(travel_angle/180*math.pi)
    # print('traveled_distance: ' + str(traveled_distance)) 
    
    q1 = dirC_to_coord (dirP_to_coord(1, pt_direction),q0,traveled_distance)
    q1 = np.array(q1)
    # print('q1: ' + str(q1)) 

    
        
    # dirP_to_coord(1, pt_direction)
    
    
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
    
    
    
    
    
    
        


