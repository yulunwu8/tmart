# This file is part of TMart.
#
# Copyright 2022 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# Surface object

import numpy as np
import pandas as pd

from scipy.interpolate import interp1d

import os.path


class SpectralSurface():
    '''Create an object to capture the spectral reflectance of surfaces when looping wavelengths. This can be used as input to reflectance in the Surface object.
    
    Arguments:

    * ``land_cover`` -- Text, currently support 'soil', 'vegetation', 'water' and 'water_chl1' ([chla]=1).
    

    Example usage::

      # Create object
      water = tmart.SpectralSurface('water_chl1')
      
      # Find spectral reflectance at 400nm
      water.wl(400)
      
      # Create a spectral surface in a numpy array 
      wl = 400 # your variable in a loop
      np.full((2, 2), water.wl(wl))
      

    '''
    def __init__(self,land_cover):
        
        file_name = os.path.join(os.path.dirname(__file__), 'ancillary', str(land_cover) + '.csv')
        
        try:
        
        
            df = pd.read_csv(file_name)
            
            wavelength = df.wl.to_numpy() * 1000
            value = df.value.to_numpy() /100
            
            self.f = interp1d(wavelength, value) # , kind='cubic')
        
        except: 
            print('Warning: failed to open ' + str(file_name))

        
    def wl(self,wavelength):
        '''
        Specify the wavelength of interest in nm and returns the reflectance of the surface. 
        '''
        
        reflectance = self.f(wavelength).item()
        if reflectance < 0: reflectance = 0
        return reflectance 





class Surface(): 
    '''Create an Surface object. 
    
    Arguments:

    * ``DEM`` -- Numpy array, Digital Elevation Model, the elevation of pixels.
    * ``reflectance`` -- Numpy array, reflectance of land or water-leaving reflectance of water, Lambertian.
    * ``isWater`` -- Numpy array, which pixels are water pixels. 1 is water, 0 is land.
    * ``cell_size`` -- The width and length of each pixel, in meters.
    * ``alignPixels`` -- Boolean. If true: the southeast corner of the pixel will take the elevation in DEM. If false: the centre of the pixel will take the elevation in DEM.

    

    Example usage::

      image_DEM = np.array([[0,0],[0,0]]) # in meters
      image_reflectance = np.array([[0.1,0.1],[0.1,0.1]]) # unitless     
      image_isWater = np.array([[1,1],[1,1]]) # 1 is water, 0 is land
      cell_size = 20_000 
      my_surface = tmart.Surface(image_DEM,image_reflectance,image_isWater,cell_size)  

    '''
    
    def __init__(self,DEM,reflectance,isWater,cell_size,alignPixels=True):
        self.DEM = DEM
        self.reflectance = reflectance
        self.isWater = isWater
        self.cell_size = cell_size
        self.alignPixels = alignPixels
        
        
        if self.DEM.shape != self.reflectance.shape:
            print('WARNING: DEM and reflectance images do not have the same shape')

        # Background parameters that can be modified using object.set_background function
        self.bg_ref = None
        self.bg_isWater = None
        self.bg_elevation = None         
        self.bg_coords = None

        # Triangulated DEM, two sets of three dimensional triangles 
        self.DEM_triangulated = None
        
        self.set_background() # test March 12, 2022
        self._triangulate_DEM()
        
        self.x_min = np.min(self.DEM_triangulated[0][:,0,:,:])
        self.x_max = np.max(self.DEM_triangulated[0][:,0,:,:])
        self.y_min = np.min(self.DEM_triangulated[0][:,1,:,:])
        self.y_max = np.max(self.DEM_triangulated[0][:,1,:,:])
        
        
        
        
    # Two Reflectance Surfaces and if water 
    def set_background(self,bg_ref=None, bg_isWater=None, bg_elevation=None, bg_coords=None):
        '''Set background information, 1 or 2 background surfaces can be set;
        If 2 surfaces: the first background is the one closer to [0,0]
        
        Arguments:

        * ``bg_ref`` -- A number or a list of two numbers. Reflectance of each background surface. Default average reflectance of the pixels.
        * ``bg_isWater`` -- A number or a list of two numbers. If each of the background surface is water. Default land.
        * ``bg_elevation`` -- A number. The elevation of the background surfaces. Default 0.
        * ``bg_coords`` -- A list of two lists of two numbers. Two XY coordinates divide the two background surfaces. Default [[0,0],[1,1]].


        Example usage::

          my_surface.set_background(bg_ref        = [0.02,0.02], # background reflectance
                                    bg_isWater    = [0,0], # if is water
                                    bg_elevation  = 0, # elevation of both background
                                    bg_coords     = [[0,0],[10,10]]) # a line dividing two background 

        '''

        # default: average reflectance, elevation 0, no water 
  
        
        # Reflectances of 2 background surfaces, first surface close to [0,0]
        if bg_ref==None: # default, average reflectance 
            self.bg_ref = [np.average(self.reflectance),np.average(self.reflectance)]  
        elif isinstance(bg_ref,list): # if list, take it 
            self.bg_ref = bg_ref
        else: # integer or float 
            self.bg_ref = [bg_ref,bg_ref]
        
        # isWater list 
        if bg_isWater==None: # default, not water 
            self.bg_isWater = [0,0]
        elif isinstance(bg_isWater,list): # if list, take it 
            self.bg_isWater = bg_isWater
        else: # integer or float 
            self.bg_isWater = [bg_isWater,bg_isWater]
        
        
        # Elevation of the background surfaces, has to be the same 
        if bg_elevation==None:
            self.bg_elevation = 0
        else:
            self.bg_elevation = bg_elevation
        
        # Two points dividing 2 background reflecting surfaces 
        if bg_coords==None: # Any value will do 
            self.bg_coords = np.array([[0,0],[1,1]])      
        else:
            self.bg_coords = np.array(bg_coords).astype(np.float32)
            # two X's can't be the same
            if self.bg_coords[0,0] == self.bg_coords[1,0]:
                self.bg_coords[1,0] = self.bg_coords[1,0] + 0.01            

        self._triangulate_DEM()


    def _triangulate_DEM(self):
        
        '''
        Make two sets of three dimensional triangles to model the surface. 
        
        Parameters
        ----------
        DEM : np array
            elevation of cells.
        cell_size : number
            width of cells.
        alignPixels: boolean
            if align with pixels (elevation will be xy positive bottom right corner of the pixel)
    
        Returns
        -------
        list
            two np arrays.
            to index the triangles: ref_tri1[point0-2, xyz:0-2, row, column]
    
        '''
        
        alignPixels=self.alignPixels
        
        if self.bg_elevation is None:
            print("====== WARNING: DEM not triangulated because of missing background elevation ======")
        else:
            
            if alignPixels: # pad two layers on left and top 
                DEM = np.pad(self.DEM,((2,1),(2,1)),'constant',constant_values = (self.bg_elevation,self.bg_elevation))
                # one elevation only  
                
            else: # pad one layer on each side 
                DEM = np.pad(self.DEM,((1,1),(1,1)),'constant',constant_values = (self.bg_elevation,self.bg_elevation))
            
            # cell_size 
            CZ = self.cell_size 
            
            ref_tri_p1 = np.array([
                # x
                np.tile(np.array(range(0,DEM.shape[1]-1)) * CZ - CZ/2,(DEM.shape[0]-1,1)),
                # y
                np.transpose([np.array(range(0,DEM.shape[0]-1)) * CZ - CZ/2]*(DEM.shape[1]-1)),
                # z
                DEM[0:DEM.shape[0]-1,0:DEM.shape[1]-1]
                ])
            
            ref_tri_p2 = np.array([
                # x
                np.tile(np.array(range(0,DEM.shape[1]-1)) * CZ - CZ/2,(DEM.shape[0]-1,1)),
                # y
                np.transpose([np.array(range(1,DEM.shape[0])) * CZ - CZ/2]*(DEM.shape[1]-1)),
                # z
                DEM[1:DEM.shape[0],0:DEM.shape[1]-1]
                ])
            
            ref_tri_p3 = np.array([
                # x
                np.tile(np.array(range(1,DEM.shape[1])) * CZ - CZ/2,(DEM.shape[0]-1,1)),
                # y
                np.transpose([np.array(range(1,DEM.shape[0])) * CZ - CZ/2]*(DEM.shape[1]-1)),
                # z
                DEM[1:DEM.shape[0],1:DEM.shape[1]]
                ])
            
            ref_tri_p4 = np.array([
                # x
                np.tile(np.array(range(1,DEM.shape[1])) * CZ - CZ/2,(DEM.shape[0]-1,1)),
                # y
                np.transpose([np.array(range(0,DEM.shape[0]-1)) * CZ - CZ/2]*(DEM.shape[1]-1)),
                # z
                DEM[0:DEM.shape[0]-1,1:DEM.shape[1]]
                ])    
            
            if alignPixels:
                ref_tri_p1[0] = ref_tri_p1[0] - CZ/2
                ref_tri_p1[1] = ref_tri_p1[1] - CZ/2
                
                ref_tri_p2[0] = ref_tri_p2[0] - CZ/2
                ref_tri_p2[1] = ref_tri_p2[1] - CZ/2
                
                ref_tri_p3[0] = ref_tri_p3[0] - CZ/2
                ref_tri_p3[1] = ref_tri_p3[1] - CZ/2
                
                ref_tri_p4[0] = ref_tri_p4[0] - CZ/2
                ref_tri_p4[1] = ref_tri_p4[1] - CZ/2
                
            ref_tri1 = np.array([ref_tri_p1, ref_tri_p2, ref_tri_p3])
            ref_tri2 = np.array([ref_tri_p1, ref_tri_p4, ref_tri_p3])
            
            self.DEM_triangulated = [ref_tri1,ref_tri2]

        












