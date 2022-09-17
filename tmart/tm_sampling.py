# This file is part of TMart.
#
# Copyright 2022 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.






# Sampling 

import random
import math
import numpy as np

from .tm_geometry import dirP_to_coord, rotation_matrix, dirC_to_dirP
from scipy.interpolate import interp1d


# Sample a direction from an isotropic distribution 
def sample_Lambertian():
    zenith = math.acos(math.sqrt(random.random())) *180/math.pi
    azimuthal = random.uniform(0,360)
    coord = dirP_to_coord(1,[zenith, azimuthal])
    
    # a list of two lists, coord & direction 
    return [coord, [zenith, azimuthal]] 

    
# Sample a scattering direction based on existing direction 
def sample_scattering(ot_mie,ot_rayleigh,pt_direction,aerosol_SPF, print_on=False): 
    
    # sum of scattering 
    ot_sum = ot_mie + ot_rayleigh

    ### Determine if mie or rayleigh 
    
    ot_random = random.uniform(0,ot_sum)
    
    # Mie 
    if ot_random <= ot_mie:
        if print_on: print ('\nMie scattering')
        type_scat = 'M'
        df_angle = aerosol_SPF.Angle.to_numpy()
        
        # sin correction 
        df_value = aerosol_SPF.Value.to_numpy() * np.sin(df_angle * math.pi/180)
        
        x_min = np.min(df_angle)
        x_max = np.max(df_angle)
        
        if x_min!=0 or x_max!=180:
            print('WARNING: Angle has to be between 0 and 180') # csv problem
 
        f2 = interp1d(df_angle, df_value, kind='cubic') 
        
        '''
        #visualize the interpolation 
        xnew = np.linspace(0, 180, num=800, endpoint=True)
        import matplotlib.pyplot as plt
        plt.plot(df_angle, df_value, 'o', xnew, f(xnew), '-', xnew, f2(xnew), '--')
        plt.legend(['data', 'linear', 'cubic'], loc='best')
        plt.show()
        '''
        
        y_max = np.max(df_value)
        y=y_max
        y_calculated=0
        
        while y>y_calculated:
            y=random.uniform(0,y_max)
            x=random.uniform(0,180)
            y_calculated = f2(x).item() # interpolate 
                
    # Rayleigh   
    else:
        if print_on: print ('\nRayleigh scattering')
        type_scat = 'R'
        
        # # Old method 
        # y_max = 0.82
        # y=y_max
        # y_calculated=0
        
        # while y>y_calculated:
        #     y=random.uniform(0,y_max)
        #     x=random.uniform(0,180)
        #     y_calculated = (3/4)*(1+(math.cos(x/180*math.pi))**2)
        #     y_calculated = y_calculated * math.sin(x/180*math.pi)
        
        # Rayleigh scattering phase function from libRadtran, line 4539 in mystic.c
        P = random.random()
        q = 8.0 * P - 4.0
        u =  (-q / 2.0 + math.sqrt (1.0 + q * q / 4.0))**(1/3)
        v = -1.0 / u
        mu = u + v
        x = math.acos(mu) * 180 / math.pi # sampled scattering angle 
        y_calculated = (3/4)*(1+(math.cos(x/180*math.pi))**2) # scattering intensity 

    # print(x) # azimuthal direction 
    sampled_direction = [x, random.uniform(0,360)]
    
    if print_on: print("Sampled_direction: " + str(sampled_direction))
    
    
    ###### Rest of the function is to rotate the sampled scattering direction towards the existing moving direction 
    
    # find an equator point using the azimuth angle of the sampled direction
    equator_point = [90, sampled_direction[1] + 90] # azimuth plus 90 to find a rotation axis 
    equator_point_C = dirP_to_coord(1,equator_point)
    
    # we rotate this equator point towards pt_direction
    # then it will be perpendicular to pt_direction 

    axis = [math.cos((pt_direction[1]+90)*math.pi/180),
            math.cos(pt_direction[1]*math.pi/180),
            0]     
    
    theta = pt_direction[0]*math.pi/180     
    equator_point_C_rotated = np.dot(rotation_matrix(axis, theta), equator_point_C)
    
    # use the new coordinates as an axis to rotate pt_direction by the scattering angle   
    coord = dirP_to_coord(1,pt_direction)
    axis2 = equator_point_C_rotated
    theta2 = sampled_direction[0]*math.pi/180 #math.pi  # zenith
    rotated = np.dot(rotation_matrix(axis2, theta2), coord)
    new_direction = dirC_to_dirP(rotated)[0:2]
    
    if print_on: print("Rotated_direction: " + str(new_direction)) 
    
    # intensity at new_direction 
    intensity = y_calculated  / math.sin(x/180*math.pi)
    
    if print_on: print ('  intensity: ' + str(intensity))
    
    return new_direction, intensity, type_scat



# Importance sampling 
def weight_impSampling(ot_mie,ot_rayleigh,angle_impSampling,aerosol_SPF, print_on=False): 
    
    # sum of scattering 
    ot_sum = ot_mie + ot_rayleigh
    ot_random = random.uniform(0,ot_sum)
    
    ### Determine if mie or rayleigh 
    
    # Mie 
    if ot_random <= ot_mie:
        if print_on: print ('\nMie scattering importance sampling')
        
        df_angle = aerosol_SPF.Angle.to_numpy()
        df_value = aerosol_SPF.Value.to_numpy() 
        
        x_min = np.min(df_angle)
        x_max = np.max(df_angle)
        
        if x_min!=0 or x_max!=180:
            print('WARNING: Angle has to be between 0 and 180') # csv problem

        # Interpolate 
        f2 = interp1d(df_angle, df_value, kind='cubic') 
        y_calculated = f2(angle_impSampling).item() 
                
    # Rayleigh   
    else:
        if print_on: print ('\nRayleigh scattering importance sampling')
        
        y_calculated = (3/4)*(1+(math.cos(angle_impSampling/180*math.pi))**2)

    # Intensity at new_direction 
    intensity = y_calculated 
    if print_on: print ('  importance sampling intensity: ' + str(intensity))

    return intensity
  
    
# Normalized to a combined SPF model 
def weight_impSampling2(ot_mie,ot_rayleigh,angle_impSampling,aerosol_SPF, print_on=False):  
    
    # sum of scattering 
    ot_sum = ot_mie + ot_rayleigh
    
    # Mie 
    df_angle = aerosol_SPF.Angle.to_numpy()
    df_value = aerosol_SPF.Value.to_numpy() 
    
    x_min = np.min(df_angle)
    x_max = np.max(df_angle)
    
    if x_min!=0 or x_max!=180:
        print('WARNING: Angle has to be between 0 and 180') # csv problem

    f2 = interp1d(df_angle, df_value, kind='cubic') 
    y_calculated_M = f2(angle_impSampling).item() # interpolate 
                
    # Rayleigh   
    y_calculated_R = (3/4)*(1+(math.cos(angle_impSampling/180*math.pi))**2)

    # intensity at new_direction 
    intensity =  y_calculated_M * (ot_mie/ot_sum) + y_calculated_R * (ot_rayleigh/ot_sum)
    
    if print_on: print ('  importance sampling intensity: ' + str(intensity))

    return intensity    
    

# # Not used anymore
# def sample_distance2scatter(ot_aerosol,ot_rayleigh,layer_height, pr):
    
#     ot_sum = ot_aerosol + ot_rayleigh
    

#     if ot_sum > 0: # if there is scattering 
#         if pr: 
#             print ("\nCalculating distance before scattered, layer_height: " + str(layer_height))
        
#         if layer_height is None:
#             print ("WARNING: atmosphere layer-height error ")
        
        
#         # attenuation coefficient due to scattering 
#         c = (ot_aerosol + ot_rayleigh )/layer_height 
        
#         R = random.random()
        
#         pt_linear_distance = - np.log(R) / c # according to Mobley light and water chapter 6
        
#         scatter_at_end = True
        
        
#     else: # if scattering coefficient is 0  
#         if pr:
#             print ('no scattering, pre-set distance')
#         pt_linear_distance = 100_000    
        
#         scatter_at_end = False
        
        
#     if pr:
#         print ('sampled_distance: ' + str(pt_linear_distance))
        
#     return pt_linear_distance, scatter_at_end




# # test sampling_scattering 


# if __name__=='__main__':
#     from Py6S.Params.atmosprofile import AtmosProfile
#     from Atmosphere import Atmosphere 
#     atm_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeSummer) # same as 6S
#     aot550 = 0.1 
#     aerosol_SPF = 'aerosol_maritime_SPF.csv'
#     aerosol_EXT = 'aerosol_maritime_EXT.csv'
#     aerosol_SSA = 'aerosol_maritime_SSA.csv'
    

#     # all wavelengths 
#     my_atm = Atmosphere(atm_profile, aot550, aerosol_SPF, aerosol_EXT, aerosol_SSA)
    
#     # 550
#     my_atm_OT, aerosol_SPF  = my_atm.wavelength(wl=550)
    
#     ot_mie = 0
#     ot_rayleigh = 1
#     pt_direction = [0,0]

#     test = sample_scattering(ot_mie,ot_rayleigh,pt_direction,aerosol_SPF)


















