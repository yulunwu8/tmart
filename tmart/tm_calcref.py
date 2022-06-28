



# Analyze results 


# When with Shadow, the individual reflectance is not reliable, total is



# Columes: pt_id, movement, type of collision, L_cox-munk, L_whitecap, L_water, L_land, L_rayleigh, L_mie, surface xyz 

# differentiate rayleigh and mie 

import numpy as np 
import pandas as pd
import sys
from copy import copy



def calc_ref(df, n_photon = None):
    '''Analyze the results of T-Mart and calculate reflectances. 
    It is possible to seperate specular and water-leaving contributions and seperate aerosol and molecule contributions. Let me know if it's useful.'
    
    Arguments:

    * ``df`` -- Results from T-Mart runs
    * ``n_photon`` -- Specify the number of photons in the run. If not specified, the number of unique pt_id will be used.
    
       * This can lead to errors when photons were fired upwards because some photons will not have pt_id.
    
    Output:

    * A list of atmospheric intrinsic reflectance, direct reflectance, environmental reflectance and total reflectance. 

    Example usage::

      R = tmart.calc_ref(results)
      print(R)

    '''
    
    print('=====================================')
    print('Calculating reflectances...')

    
    
    dfpd = pd.DataFrame(data=df, columns = ['pt_id', 'movement', 'type_collision', 'L_coxmunk', 
                                            'L_whitecap', 'L_water', 'L_land', 'L_rayleigh', 'L_mie','x','y','z','shadowed'])
       
    
    dfpd.pt_id =        pd.to_numeric( dfpd.pt_id)
    dfpd.movement =     pd.to_numeric( dfpd.movement)
    dfpd.L_coxmunk =    pd.to_numeric( dfpd.L_coxmunk)
    dfpd.L_whitecap =   pd.to_numeric( dfpd.L_whitecap)
    dfpd.L_water =      pd.to_numeric( dfpd.L_water)
    dfpd.L_land =       pd.to_numeric( dfpd.L_land)
    dfpd.L_rayleigh =   pd.to_numeric( dfpd.L_rayleigh)
    dfpd.L_mie =        pd.to_numeric( dfpd.L_mie)
    dfpd.x =            pd.to_numeric( dfpd.x)
    dfpd.y =            pd.to_numeric( dfpd.y)
    dfpd.z =            pd.to_numeric( dfpd.z)
    dfpd.shadowed =     pd.to_numeric( dfpd.shadowed)
    
    dfpd['env'] = 0 # 1 is environment reflectance 
    
    unique_pt_id = (np.unique(dfpd.pt_id))
    
    if n_photon == None:
        n_photon = unique_pt_id.shape[0]

    
    dfpd.dtypes
    
    ###
    
    # number of unique pt_id
    pt_id_counts = dfpd.pt_id.value_counts()
    
    # isolate those that only have one
    pt_ids_1 = pt_id_counts.index[pt_id_counts == 1].tolist()
    
    # df_output = dfpd.iloc[pt_ids_1] 
    
    df_output = dfpd[dfpd['pt_id'].isin(pt_ids_1)] # append to this 
    
    
    
    
    pt_ids_2 = pt_id_counts.index[pt_id_counts > 1].tolist() # pt_ids more than 1 rows 
    
    
    
    pt_ids = pt_ids_2
    
    # pt_ids.sort()
    
    
    
    
    # slower method 
    
    # df_output = pd.DataFrame()
    # pt_ids = pd.unique(dfpd.pt_id)
    
    
    for pt_id in pt_ids:
        # print(pt_id)
    
        # for every pt_id, we only need to process the ones with more than 1 row!!! 
        # with 1 rows copy to new data frame 
        
        # pt_id = 196
    
        # data frame with this pt_id 
        pt = dfpd[dfpd.pt_id  == pt_id]
    
    
        moves = pt.movement # loop through this
        
        if not moves.is_monotonic_increasing: # check if sorted 
            sys.exit('pt movement has to be sorted')
        
        
        # Find the move with W, L or Ws collision, then add all after to them and stop the loop
        
        for move in moves:
    
            # move = 0
            
            pt_movement = copy(pt[pt.movement == move]) # entire row 
            
            t_c = pt_movement.type_collision.values[0] # type of collision 
            type(t_c)
            
            if t_c=='W' or t_c=='L':
                # print('Adding all after to L_whitecap, L_water, L_land')
                
                
                if sum( pt.iloc[move,4:7]) == 0: # if surface is black 
                    break
                
                # pt.iloc[move+1:,3:9] # all below
                
                ### Old way before introducing shadow
                # sum_after = np.sum(pt.iloc[move+1:,3:9].values)
                
                pt_mov_after = pt.iloc[move+1:,:]
                pt_mov_after_nonShadow = pt_mov_after[pt_mov_after.shadowed==0]
                sum_after = np.sum(pt_mov_after_nonShadow.iloc[:,3:9].values)
                
                
                # total of the single row
                total = pt_movement.L_whitecap.values[0] + pt_movement.L_water.values[0] + pt_movement.L_land.values[0]
                
                # ratio 
                r_wc = pt_movement.L_whitecap.values[0] / total 
                r_water = pt_movement.L_water.values[0] / total
                r_land = pt_movement.L_land.values[0] / total 
                
                # We calculate ratio before setting 'total' to 0
                if pt_movement.shadowed.values[0] == 1: total = 0
                
                total_new = total + sum_after
            
                pt_movement.L_whitecap = total_new * r_wc
                pt_movement.L_water = total_new * r_water
                pt_movement.L_land = total_new * r_land
                
                if move>0: pt_movement.env = 1
                
                # df_output = df_output.append(pt_movement)
                
                df_output = pd.concat([df_output,pt_movement])
                
                
                break 
                
            if t_c=='Ws':    
                # print('Adding all after to L_coxmunk')
                
                ### Old way
                # sum_after = np.sum(pt.iloc[move+1:,3:9].values)    
                
                pt_mov_after = pt.iloc[move+1:,:]
                pt_mov_after_nonShadow = pt_mov_after[pt_mov_after.shadowed==0]
                sum_after = np.sum(pt_mov_after_nonShadow.iloc[:,3:9].values)          
                
                if pt_movement.shadowed.values[0] == 1: 
                    total = 0
                else:
                    total = pt_movement.L_coxmunk.values[0]
                
                pt_movement.L_coxmunk = total + sum_after
                
                if move>0: pt_movement.env = 1
                
                # df_output = df_output.append(pt_movement)
                
                df_output = pd.concat([df_output,pt_movement])
                
                break 
            
            # df_output = df_output.append(pt_movement)
            
            # If rayleigh or mie 
            df_output = pd.concat([df_output,pt_movement])
        
    
    
    # End of loop 
    
    
    # if mov == 0 and type_collision == W or L, then all after add to L_whitecap, L_water, L_land
    # if type_collision == Ws, to L_coxmunk instead 
    
    # if mov > 0, do the same except env = 1 
    # replace the used ones with 0!!!
    
    
    
    
    R_atm = np.sum(df_output.iloc[:,7:9].values) / n_photon
    
    R_dir = np.sum(df_output[df_output.env==0].iloc[:,3:7].values) / n_photon
    
    R_env = np.sum(df_output[df_output.env==1].iloc[:,3:7].values) / n_photon
    
    R_total = R_atm + R_dir + R_env
    
    
    R_output = {'R_atm':R_atm,
                'R_dir':R_dir,
                'R_env':R_env,
                'R_total':R_total}
    
    
    return R_output
    
     










