# This file is part of TMart.
#
# Copyright 2022 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# Geometry 


import math
import numpy as np
import sys


def rad(degree): 
    return degree/180*math.pi



def dirP_to_coord(distance, direction):
    '''
    convert polar coordinates and distance to coordinates 
    
    
    Parameters
    ----------
    distance : integer or float 
        Vector: the linear distance or speed.
    direction : list of two integer or float 
        First zenith: up positive Z 0, positive Y 90. 0-180
        Second azimuthal: positive X 0, positive Y 90. 0-360

    Returns
    -------
    coords : list
        Vector in 3 dimensions: X, Y, and Z.

    '''
    # Two inputs: vector

    coords = np.array([0,0,0],float)  
    #math.sqrt(coords[0]**2 + coords[1]**2 +coords[2]**2) 
    
    coords[0] = distance * math.sin(rad(direction[0])) * math.cos(rad(direction[1])) # x
    coords[1] = distance * math.sin(rad(direction[0])) * math.sin(rad(direction[1])) # y
    coords[2] = distance * math.cos(rad(direction[0])) #z
    return coords



def angle_3d(a,b,c):
    '''
    
    calculate the angle between 3 points in 3d 
    source: https://stackoverflow.com/questions/35176451/python-code-to-calculate-angle-between-three-point-using-their-3d-coordinates
    
    '''
    
    a,b,c = np.array(a),np.array(b),np.array(c)

    ba = a - b
    bc = c - b
    
    
    if np.linalg.norm(ba)==0:
        print(a)
        print(b)
        print(c)
        sys.exit("error1 angle_3d")
        
    if np.linalg.norm(bc)==0:
        print(a)
        print(b)
        print(c)
        sys.exit("error2 angle_3d")        
    
    
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    
    # np.rad2deg(2.1034) # useful numpy operations 
    # np.deg2rad(120.5159)
    
    return np.degrees(angle)



def linear_distance(pt1, pt2):
    '''
    each point should be a list of 3 coordinates: X, Y, Z
 
    '''
    distance = ((pt2[0]-pt1[0])**2 + (pt2[1]-pt1[1])**2 + (pt2[2]-pt1[2])**2)**(1/2)
    return distance


def dirC_to_dirP(xyz): 
    #takes list xyz (single coord), return theta, phi, r 
    x       = xyz[0]
    y       = xyz[1]
    z       = xyz[2]
    r       =  (x*x + y*y + z*z)**0.5
    theta   =  math.acos(z/r)*180/ math.pi #to degrees
    phi     =  math.atan2(y,x)*180/ math.pi
    
    if phi <0:
        phi = phi +360
    
    return [theta,phi,r]



def rotation_matrix(axis, theta): 
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    source: https://stackoverflow.com/questions/6802577/rotation-of-3d-vector
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])




def dirC_to_coord (direction,q0,linear_length): 
    '''
    scale a cartesian direction to coordinates with a starting point and a length 

    '''
    
    linear_distance = ((direction[0]-q0[0])**2 + 
                       (direction[1]-q0[1])**2 + 
                       (direction[2]-q0[2])**2)**(1/2)
    
    
    linear_distance = ((direction[0])**2 + 
                       (direction[1])**2 + 
                       (direction[2])**2)**(1/2)
    
    scale = linear_length / linear_distance
    
    scaled = [q0[0]+(direction[0]) * scale,
              q0[1]+(direction[1]) * scale,
              q0[2]+(direction[2]) * scale]
    
    return scaled






