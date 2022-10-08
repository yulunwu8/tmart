# This file is part of TMart.
#
# Copyright 2022 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.







# Analyze the output of TMart and differentiate direct, env and atm intrinsic reflectances  

import numpy as np 
import pandas as pd
import sys
from copy import copy

def calc_ref(df, n_photon = None, detail = False, total_only = False):
    '''Analyze the results of T-Mart and calculate reflectances. 
    
    Arguments:

    * ``df`` -- Results from T-Mart runs
    * ``n_photon`` -- Specify the number of photons in the run when firing the photon upwards. 
       * If not specified, the number of unique pt_id will be used. This can lead to errors when photons were fired upwards because some photons will not have pt_id.
    * ``detail`` -- differentiate coxmunk, whitecap, water-leaving and land contributions
    * ``total_only`` -- only return total TOA reflectance, this is much faster but gives less details. 
    
    Output:

    * A list of atmospheric intrinsic reflectance, direct reflectance, environmental reflectance and total reflectance. 

    Example usage::

      R = tmart.calc_ref(results)
      for k, v in R.items():
          print(k, '\t ' , v)

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
    
    dfpd['env'] = 0 # mark environment reflectance as 1 
    
    unique_pt_id = (np.unique(dfpd.pt_id))
    
    if n_photon == None:
        n_photon = unique_pt_id.shape[0]
    
    
    ### Total reflectance only 
    
    if total_only:
        R_total = np.sum(dfpd.L_coxmunk) + np.sum(dfpd.L_whitecap) + np.sum(dfpd.L_water) + np.sum(dfpd.L_land) + np.sum(dfpd.L_rayleigh) + np.sum(dfpd.L_mie) 
        R_total = R_total / n_photon
        return R_total
    

    
    ### Differentiate atm, dir and env reflectances 
    
    # Number of unique pt_id
    pt_id_counts = dfpd.pt_id.value_counts()   
    
    # Isolate those that only have one pt_id to speed up calculation 
    pt_ids_1 = pt_id_counts.index[pt_id_counts == 1].tolist()
    
    # Main output data frame 
    df_output = dfpd[dfpd['pt_id'].isin(pt_ids_1)] # append to this 
    
    # Extract the ones with more than one pt_id, meaning multiple scattering events 
    pt_ids_2 = pt_id_counts.index[pt_id_counts > 1].tolist() # pt_ids more than 1 rows 
    
    '''testing 
    
    pt_id = 1
    
    '''
    
    # Loop through pt_ids, differentiate reflectances 
    for pt_id in pt_ids_2:
    
        # Data frame with this pt_id 
        pt = dfpd[dfpd.pt_id  == pt_id]
        
        # Each movement of the photon 
        moves = pt.movement 
        
        if not moves.is_monotonic_increasing: # check if sorted 
            sys.exit('pt movement has to be sorted')
        
        # Find the move with W, L or Ws collision, then add all after to them and stop the loop
        for move in moves:
            
            # Entire row 
            pt_movement = copy(pt[pt.movement == move]) 
            
            # Type of collision 
            t_c = pt_movement.type_collision.values[0] 
            
            if t_c=='W' or t_c=='L':
                # Adding all after to L_whitecap, L_water, L_land
                
                # If surface is black, quit loop to speed up calculation 
                if sum( pt.iloc[move,4:7]) == 0: 
                    break
                
                # Identify all movements after that contribute to this 
                pt_mov_after = pt.iloc[move+1:,:]
                pt_mov_after_nonShadow = pt_mov_after[pt_mov_after.shadowed==0]
                sum_after = np.sum(pt_mov_after_nonShadow.iloc[:,3:9].values)
                
                # Total of the single row
                total = pt_movement.L_whitecap.values[0] + pt_movement.L_water.values[0] + pt_movement.L_land.values[0]
                
                # Ratios 
                r_wc = pt_movement.L_whitecap.values[0] / total 
                r_water = pt_movement.L_water.values[0] / total
                r_land = pt_movement.L_land.values[0] / total 
                
                # We calculate ratio before setting shadowed 'total' to 0
                if pt_movement.shadowed.values[0] == 1: total = 0
                
                total_new = total + sum_after
            
                pt_movement.L_whitecap = total_new * r_wc
                pt_movement.L_water = total_new * r_water
                pt_movement.L_land = total_new * r_land
                
                # move>0 means at least one atmospheric scattering happened
                if move>0: pt_movement.env = 1
                
                # Add to the main data frame 
                df_output = pd.concat([df_output,pt_movement])
                
                break 
                
            if t_c=='Ws':    
                # Adding all after to L_coxmunk
                
                pt_mov_after = pt.iloc[move+1:,:]
                pt_mov_after_nonShadow = pt_mov_after[pt_mov_after.shadowed==0]
                sum_after = np.sum(pt_mov_after_nonShadow.iloc[:,3:9].values)          
                
                if pt_movement.shadowed.values[0] == 1: 
                    total = 0
                else:
                    total = pt_movement.L_coxmunk.values[0]
                
                pt_movement.L_coxmunk = total + sum_after
                
                if move>0: pt_movement.env = 1
                
                df_output = pd.concat([df_output,pt_movement])
                
                break 
            
            # If rayleigh or mie 
            df_output = pd.concat([df_output,pt_movement])
    
    
    R_atm = np.sum(df_output.iloc[:,7:9].values) / n_photon
    R_dir = np.sum(df_output[df_output.env==0].iloc[:,3:7].values) / n_photon
    R_env = np.sum(df_output[df_output.env==1].iloc[:,3:7].values) / n_photon
    R_total = R_atm + R_dir + R_env
    
    
    if detail:
    
        R_dir_coxmunk = np.sum(df_output[df_output.env==0].iloc[:,3].values) / n_photon
        R_dir_whitecap = np.sum(df_output[df_output.env==0].iloc[:,4].values) / n_photon
        R_dir_water = np.sum(df_output[df_output.env==0].iloc[:,5].values) / n_photon
        R_dir_land = np.sum(df_output[df_output.env==0].iloc[:,6].values) / n_photon
        
        R_env_coxmunk = np.sum(df_output[df_output.env==1].iloc[:,3].values) / n_photon
        R_env_whitecap = np.sum(df_output[df_output.env==1].iloc[:,4].values) / n_photon
        R_env_water = np.sum(df_output[df_output.env==1].iloc[:,5].values) / n_photon
        R_env_land = np.sum(df_output[df_output.env==1].iloc[:,6].values) / n_photon        
        
        R_output = {'R_atm':R_atm,
                    'R_dir':R_dir,
                    '_R_dir_coxmunk':R_dir_coxmunk,
                    '_R_dir_whitecap':R_dir_whitecap,
                    '_R_dir_water':R_dir_water,
                    '_R_dir_land':R_dir_land,
                    'R_env':R_env,
                    '_R_env_coxmunk':R_env_coxmunk,
                    '_R_env_whitecap':R_env_whitecap,
                    '_R_env_water':R_env_water,
                    '_R_env_land':R_env_land,
                    'R_total':R_total}
    
    else:      
        R_output = {'R_atm':R_atm,
                    'R_dir':R_dir,
                    'R_env':R_env,
                    'R_total':R_total}
    
    
    return R_output
    
     








