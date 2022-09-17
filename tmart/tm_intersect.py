# This file is part of TMart.
#
# Copyright 2022 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# intersect


import pandas as pd
import numpy as np
import sys

from .tm_geometry import angle_3d, linear_distance
from copy import copy



# new find_atm two scattering OTs for travelling through multiple layers 
def find_atm2(atm_profile,q1):
    
    z = q1[2]/1000
    
    if z < 0:
        print("ERROR: z < 0")
        return None
    
    # print(atm_profile)
    
    alt_top_minus_z = atm_profile[:,1] - z
    # print('alt_top_minus_z: ' + str(alt_top_minus_z))
    
    
    # index of minimum positive difference 
    idx = np.where(alt_top_minus_z >= 0, alt_top_minus_z, np.inf).argmin()
    # print('Index: ' + str(idx))
    
    
    # print(atm_profile[idx,:])
    
    ot_rayleigh = atm_profile[idx,3]
    ot_mie = atm_profile[idx,4]
    
    return ot_rayleigh, ot_mie
    




# Testing intersecting triangles 
# 2: test boxes first 
# This one is used in the code because it is a lot faster 

def intersect_line_DEMtri2(q0, q1, DEM_tri, print_on = False):
    '''

    Parameters
    ----------
    q0 : list
        starting point of the line.
    q1 : list
        ending point.
    DEM_tri : list
        a list of numpy arrays.

    Returns
    -------
    intersect_tri : pandas dataframe
        collision coordinates, direction normal to the surface
        and linear distance to the starting point.

    '''
    
    # Manual switch 
    print_on = False 
    
    
    ### Identify boxes intersecting the line first     

    tri_x = DEM_tri[0][:,0] 
    tri_y = DEM_tri[0][:,1] 
    tri_z = DEM_tri[0][:,2]
    
    
    tri_x_min = np.min(tri_x,axis=0)
    tri_x_max = np.max(tri_x,axis=0)

    tri_y_min = np.min(tri_y,axis=0)  
    tri_y_max = np.max(tri_y,axis=0)    
    
    
    # Find the box within max z and min z, project it onto XY and select the pixels 
       
    pt_direction_c = [q1[0]-q0[0],q1[1]-q0[1],q1[2]-q0[2]]
    
    z_max = np.max(tri_z) + 0.1 
    z_min = np.min(tri_z) - 0.1
    
    pt_direction_c = q1 - q0
    
    # make box 
    q_z_max = q0 + ( (z_max - q0[2]) / pt_direction_c[2] ) * pt_direction_c
    q_z_min = q0 + ( (z_min - q0[2]) / pt_direction_c[2] ) * pt_direction_c
    
    # First treat all triangles as boxes (max/min XYZ)!!!
    
    box_x_min = min(q_z_max[0],q_z_min[0])
    box_x_max = max(q_z_max[0],q_z_min[0])
    box_y_min = min(q_z_max[1],q_z_min[1])
    box_y_max = max(q_z_max[1],q_z_min[1])
    
    # boxes that cross z
    crossing_z = np.logical_and(np.logical_and(tri_x_min < box_x_max, tri_x_max > box_x_min),
                                np.logical_and(tri_y_min < box_y_max, tri_y_max > box_y_min))
    
    # crossing is for all xyz 
    
    if np.sum(crossing_z) <= 1:
        crossing = crossing_z
        
    else:
    
        if q0[0] == q1[0] and q0[1] == q1[1]:
            if print_on: print('Vertical move')
            crossing = crossing_z
            
        elif q0[0] == q1[0]:
            if print_on: print('Equal x, move along y')
            crossing_y = np.logical_and(tri_x_min < q0[0], tri_x_max > q0[0])
            crossing = np.logical_and(crossing_z, crossing_y)
            
            
        elif q0[1] == q1[1]:
            if print_on: print('Equal y, move along x')    
            crossing_x = np.logical_and(tri_y_min < q0[1], tri_y_max > q0[1])
            crossing = np.logical_and(crossing_z, crossing_x)
            
            
        else:
            if print_on: print('Move along both x and y')
            
            # test if intersects the top or bottom of the square
            
            # equation on the x-y plane
            slope = (q1[1] - q0[1]) / (q1[0] - q0[0])
            y_itcp = q1[1] - q1[0]*slope
            
            
            # crossing the upper boundary (more negative y)
            x_itcp_y_min = (tri_y_min - y_itcp) / slope
            
            # crossing the lower boundary (more positive y)
            x_itcp_y_max = (tri_y_max - y_itcp) / slope
            
            y_itcp_x_min = tri_x_min * slope + y_itcp
            y_itcp_x_max = tri_x_max * slope + y_itcp

            
            crossing_upper = np.logical_and(tri_x_min < x_itcp_y_min, x_itcp_y_min < tri_x_max)
            crossing_lower = np.logical_and(tri_x_min < x_itcp_y_max, x_itcp_y_max < tri_x_max)
            crossing_x = np.logical_or(crossing_upper,crossing_lower) 
            
            
            crossing_left = np.logical_and(tri_y_min < y_itcp_x_min, y_itcp_x_min < tri_y_max)
            crossing_right = np.logical_and(tri_y_min < y_itcp_x_max, y_itcp_x_min < y_itcp_x_max)
            crossing_y = np.logical_or(crossing_left,crossing_right) 
            
            
            crossing_xy = np.logical_or(crossing_x,crossing_y) 
            crossing = np.logical_and(crossing_z, crossing_xy)
            
            
            # np.where(crossing_z)
            # np.where(crossing_xy)
            
            

    # intersecting triangles 
    intersect_tri = pd.DataFrame() 
    
    for tri in DEM_tri: # tri is a set of triangle, ref_tri1[point0-2, xyz:0-2, row, column]
    
        # tri[:,:,crossing]
        for i in range(0,tri[:,:,crossing].shape[2]):

            # print (i)
            
            p0 = tri[:,:,crossing] [0,:,i]
            p1 = tri[:,:,crossing] [1,:,i]
            p2 = tri[:,:,crossing] [2,:,i]
                
            # Extract coordinates of intersection if exists 
            intersect = _intersect_line_triangle(q0,q1,p0,p1,p2)

            
            if (intersect is not None):             
                
                if print_on:
                    print("Line intersecting triangle at: " + str(intersect))
                
                N = np.cross(p1-p0, p2-p0) 
                
                # test normal if in the same direction as the incoming line 
                # if not, reverse normal direction 
                if angle_3d(q0,intersect,(intersect + N)) > 90: N = -N
                
                intersect_tri_temp = pd.DataFrame({
                    # coordinates and direction 
                    'X': [intersect[0]],
                    'Y': [intersect[1]],
                    'Z': [intersect[2]],
                    'N_X': [N[0]],
                    'N_Y': [N[1]],
                    'N_Z': [N[2]]
                    })
                
                intersect_tri_temp['linear_distance'] = linear_distance(
                    intersect_tri_temp.iloc[0,0:3].tolist(), 
                    q0)
                
                # intersect_tri = intersect_tri.append(intersect_tri_temp) # append to concat on Feb 4, 2022
                intersect_tri = pd.concat([intersect_tri,intersect_tri_temp])
      
    intersect_tri.reset_index(drop=True, inplace=True)
    return intersect_tri





def _intersect_line_triangle(q1,q2,p1,p2,p3):
    '''

    Parameters
    ----------
    q1 : list
        Coordinates of one end of the line: X, Y, Z.
    q2 : list
        Coordinates of the other end of the line.
    p1 : list
        Coordinates of the first angle of the triangle.
    p2 : list
        Coordinates of the second angle of the triangle.
    p3 : list
        Coordinates of the third angle of the triangle.

    Returns
    -------
    Coordinates of the intersection if they intersect, or 
    None if they don't.


    '''
    # from https://stackoverflow.com/questions/42740765/intersection-between-line-and-triangle-in-3d
    def signed_tetra_volume(a,b,c,d):
        return np.sign(np.dot(np.cross(b-a,c-a),d-a)/6.0)

    s1 = signed_tetra_volume(q1,p1,p2,p3)
    s2 = signed_tetra_volume(q2,p1,p2,p3)

    if s1 != s2:
        s3 = signed_tetra_volume(q1,q2,p1,p2)
        s4 = signed_tetra_volume(q1,q2,p2,p3)
        s5 = signed_tetra_volume(q1,q2,p3,p1)
        if s3 == s4 and s4 == s5:
            n = np.cross(p2-p1,p3-p1)
            t = np.dot(p1-q1,n) / np.dot(q2-q1,n)
            return q1 + t * (q2-q1)
    return None





def reflectance_intersect(q_collision, image_reflectance, cell_size, bg_ref, bg_coords):
    '''
    Find the reflectance of a triangle when collision happens
    June 10, 2021: also if is water 

    Parameters
    ----------
    q_collision : TYPE
        DESCRIPTION.
    image_reflectance : TYPE
        DESCRIPTION.
    cell_size : TYPE
        DESCRIPTION.
    b1_ref : TYPE
        background reflectance 1.
    b2_ref : TYPE
        background reflectance 2.
    b1_a : TYPE
        first point dividing b1_ref and b2_ref.
    b1_b : TYPE
        2nd point.

    Returns
    -------
    triangle_reflectance : TYPE
        DESCRIPTION.

    '''
    
    b1_ref = bg_ref[0]
    b2_ref = bg_ref[1]
    
    image_shape = image_reflectance.shape

    b1_a, b1_b = np.array(bg_coords[0]), np.array(bg_coords[1])
    
    triangle_y = int(q_collision[1]//cell_size)
    triangle_x = int(q_collision[0]//cell_size)
         
    line = b1_a - b1_b
    slope = line[1]/line[0]  # y = slope * x + x_intercept 
    x_intercept = b1_a[1] - slope *b1_a[0]   
    
    # if on the padding triangles 
    
    if (triangle_y < 0 or triangle_x < 0 or 
        triangle_y >= image_shape[0] or triangle_x >= image_shape[1]):
    
        above_line = q_collision[1] >= (q_collision[0] * slope + x_intercept) 
        
        if x_intercept>=0:
            #print('below the line is b1')
            
            if above_line:
                triangle_reflectance = b2_ref
            else:
                triangle_reflectance = b1_ref
                
        else:
            #print('above the line is b1')
            
            if above_line:
                triangle_reflectance = b1_ref
            else:
                triangle_reflectance = b2_ref
        
    else: # if within the original image 
    
        triangle_reflectance =  image_reflectance[triangle_y, triangle_x]
    
    return triangle_reflectance



def reflectance_background(q_collision,bg_ref, bg_coords):
    # Find reflectance when outside the triangles 
    # June 10, 2021: also if is water 
    
    b1_ref = bg_ref[0]
    b2_ref = bg_ref[1]
    
    b1_a, b1_b = np.array(bg_coords[0]), np.array(bg_coords[1])
    line = b1_a - b1_b
    slope = line[1]/line[0]  # y = slope * x + x_intercept 
    x_intercept = b1_a[1] - slope *b1_a[0]  
    
    above_line = q_collision[1] >= (q_collision[0] * slope + x_intercept) 
        
    if x_intercept>=0:
        #print('below the line is b1')
        
        if above_line:
            triangle_reflectance = b2_ref
        else:
            triangle_reflectance = b1_ref
            
    else:
        #print('above the line is b1')
        
        if above_line:
            triangle_reflectance = b1_ref
        else:
            triangle_reflectance = b2_ref
    return triangle_reflectance


def intersect_background(q0,q1,bg_elevation):
    '''
    Find the XYZ of q0-q1 at Z 

    '''

    ratio = (q0[2]-bg_elevation) / (q0[2]-q1[2])
    
    x = q0[0] - ((q0[0]-q1[0])) * ratio
    y = q0[1] - ((q0[1]-q1[1])) * ratio

    return [x,y,bg_elevation]


# test if line intersects atm boundaries 
def intersect_line_boundary(q0,q1,atm_profile,print_on=False):
    '''
    

    return the closest layer 


    '''
    
    '''
    ### testing 
    q0 = [0,0,110_000]
    q1 = [0,0,40_000]
    q0 = [0,0,90_000]
    q1 = [0,0,110_000]    
    ###
    '''
    
    z0, z1 = q0[2]/1000,q1[2]/1000
    
    if z0==z1:
        sys.exit("intersect_line_boundary: z0==z1")
    
    going_up = z1>z0
    # print(going_up)
    
    cross_b = False # if cross a boundary 
    out = False # if out
    
    # output 
    intersect_b = pd.DataFrame() # default output if no crossing boundary 
    
    
    if going_up:
        diff = atm_profile.Alt_bottom - z0
        
        # print('======')
        # print('z0: ' + str(z0))
        # print('atm_profile.Alt_bottom: ' + str(atm_profile.Alt_bottom))
        # print('diff: ' + str(diff))

        # print('diff>0: ' + str(diff>0))
        
        # print('any: ' + str(any(diff>0)))
        
        if any(diff>0): # at least one alt_bottom above z0
            layer_idx = diff[diff>0].idxmin()
            
            # print(layer_idx)
            
            layer_bottom = atm_profile.Alt_bottom[layer_idx]
            # print(layer_bottom)
            
            if (z1 >= layer_bottom):
                # append info 
                cross_b = True 
                altitude = layer_bottom # altitude of crossing boundary 
        else: # all alt_bottom below z0 
            if (z1>=100 and going_up):
                out = True 
                cross_b = True
                altitude = 100
                layer_idx = atm_profile.Alt_bottom.idxmax() # return top layer because OTs are not used 
            
    
    
    else: # going down 
        diff = z0 - atm_profile.Alt_top
        # print('======')
        # print('z0: ' + str(z0))
        # print('atm_profile.Alt_top: ' + str(atm_profile.Alt_top))
        # print('diff: ' + str(diff))

        # print('diff>0: ' + str(diff>0))
        
        # print('any: ' + str(any(diff>0)))
        
        if any(diff>0): # no alt_top above q0 --> no need to edit atm 

            layer_idx = diff[diff>0].idxmin()
            
            # print(layer_idx)
            
            layer_top = atm_profile.Alt_top[layer_idx]
            # print(layer_top)
            
            if(z1 <= layer_top):
                # append info 
                cross_b = True 
                altitude = layer_top # altitude of crossing boundary 
            
    if cross_b: 
        intersect_b = copy(atm_profile.iloc[layer_idx])
        intersect_b['out'] = False
        
        # where crossing boundary  
        xyz = intersect_background(q0,q1,altitude*1000) + [altitude*1000]
        
        intersect_b['X'] = xyz[0]
        intersect_b['Y'] = xyz[1]
        intersect_b['Z'] = xyz[2]
        
        intersect_b['linear_distance'] = linear_distance(q0,xyz)

        
        if out:
            intersect_b['out'] = True
    
    return intersect_b 







'''
# Not used anymore, travel layer by layer 
def find_atm(q0, pt_direction, atm_profile_wl, print_on = False): 
    # Find a pohton in the atmosphere, retrieve OTs 
    
    # Test if the photon alt is equal to any of the upper boundaries 
    equal_top = np.array( q0[2]/1000 == atm_profile_wl.Alt_top )
      
    # At boundary 
    if sum(equal_top): 
    
        atm_index = np.where(equal_top )[0].astype(int)[0]
        
        # Look at moving direction, if moving up, choose the layer above 
        if pt_direction[0]<90:  
            atm_index = atm_index + 1
        elif pt_direction[0]==90: 
            print("\n====== WARNING: photon cannot travel horizontally at a boundary ======")  
            atm_index = None
        
        if print_on:
            print('\nPhoton alt at a boundary, atm_index: ' + str(atm_index))
        
        ot_abs = float(atm_profile_wl.iloc[atm_index].ot_abs)
        ot_rayleigh = float(atm_profile_wl.iloc[atm_index].ot_rayleigh)
        ot_mie = float(atm_profile_wl.iloc[atm_index].ot_mie)        
        
        
    # Not at boundary   
    else:
        
        above_bottom = np.array( q0[2]/1000 > atm_profile_wl.Alt_bottom )
        below_top = np.array( q0[2]/1000 < atm_profile_wl.Alt_top )    
        
        above_or_below = above_bottom & below_top
        
        # if print_on:
        #     print('\n Above_bottom: ')
        #     print(above_bottom)
        #     print('\n Below_top: ')
        #     print(below_top)
        #     print('\n Above_or_below: ')
        #     print(above_or_below)

        # Within a layer 
        if sum(above_or_below): 
            
            atm_index = np.where(above_bottom & below_top)[0].astype(int)[0]
            
            if print_on:
                print('\nPhoton alt within a layer, atm_index: ' + str(atm_index))
            
            ot_abs = float(atm_profile_wl.iloc[atm_index].ot_abs)
            ot_rayleigh = float(atm_profile_wl.iloc[atm_index].ot_rayleigh)
            ot_mie = float(atm_profile_wl.iloc[atm_index].ot_mie)

        # Not at boundary nor within a layer 
        else: 
            # Above all layers 
            if q0[2]/1000 > np.max(atm_profile_wl.Alt_top):
                
                if print_on:
                    print ("\nPhoton alt is above max, setting OTs as 0 ")
                    
                ot_abs, ot_rayleigh, ot_mie = 0, 0, 0
            
            # Unknown 
            else:
                print("====== WARNING: error locating photon in the atm ======")  
                
    if print_on:
        print("ot_abs: {:.2e}, ot_rayleigh: {:.2e}, ot_mie: {:.2e}".format(ot_abs, ot_rayleigh, ot_mie))
    
    return ot_abs, ot_rayleigh, ot_mie
'''





'''
# Testing intersecting triangles, this function is not used 
def intersect_line_DEMtri(q0, q1, DEM_tri, print_on = False): # not used 

    # Parameters
    # ----------
    # q0 : list
    #     starting point of the line.
    # q1 : list
    #     ending point.
    # DEM_tri : list
    #     a list of numpy arrays.

    # Returns
    # -------
    # intersect_tri : pandas dataframe
    #     collision coordinates, direction normal to the surface
    #     and linear distance to the starting point.
    
    # intersecting triangles 
    intersect_tri = pd.DataFrame() 
    
    for tri in DEM_tri: 
    
        for row in range (0, tri.shape[2]):
            for col in range (0, tri.shape[3]):
                #print (row, col)
        
                p0 = tri[0,:,row,col] 
                p1 = tri[1,:,row,col]
                p2 = tri[2,:,row,col]
                
                # Extract coordinates of intersection if exists 
                intersect = _intersect_line_triangle(q0,q1,p0,p1,p2)
                if (intersect is not None): 
                    
                    
                    if print_on:
                        print("Line intersecting triangle at: " + str(intersect))
                    
                    N = np.cross(p1-p0, p2-p0) #maybe faster than the line above 
                    

                    if angle_3d(q0,intersect,N) > 90: N = -N
                    
                    intersect_tri_temp = pd.DataFrame({

                        'X': [intersect[0]],
                        'Y': [intersect[1]],
                        'Z': [intersect[2]],
                        'N_X': [N[0]],
                        'N_Y': [N[1]],
                        'N_Z': [N[2]]
                        })
                    
                    intersect_tri_temp['linear_distance'] = linear_distance(
                        intersect_tri_temp.iloc[0,0:3].tolist(), 
                        q0)
                    
                    # intersect_tri = intersect_tri.append(intersect_tri_temp) # append to concat on Feb 4, 2022
                    intersect_tri = pd.concat([intersect_tri,intersect_tri_temp])
      
    intersect_tri.reset_index(drop=True, inplace=True)
    return intersect_tri
'''









