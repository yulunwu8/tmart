# Surface object

import numpy as np
import pandas as pd

from scipy.interpolate import interp1d




class spectral_surface():
    '''Create an object to capture the spectral reflectance of surfaces when looping wavelengths.
    
    Arguments:

    * ``land_cover`` -- Text, currently support soil (general), vegetation (general), water (general), water_chl1 ([chla]=1).
    

    Example usage::

      # Create object
      water = tmart.spectral_surface('water_chl1')
      
      # Find spectral reflectance at 400nm
      water.wl(400)

    '''
    def __init__(self,land_cover):
        file_name = 'tmart/ancillary/' + str(land_cover) + '.csv'
        
        try:
        
        
            df = pd.read_csv(file_name)
            
            wavelength = df.wl.to_numpy() * 1000
            value = df.value.to_numpy() /100
            
            self.f = interp1d(wavelength, value, kind='cubic')
        
        except: 
            print('Warning: failed to open ' + str(file_name))

        
    def wl(self,wavelength):
        '''
        See above.
        '''
        
        reflectance = self.f(wavelength).item()
        if reflectance < 0: reflectance = 0
        return reflectance 









class Surface(): 
    '''Create an Atmosphere object. 
    
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
    
    # add water, wavelength, function 
    # add DEM/reflectance not the same coordinates - Use DEM, find reflectance from image
    # add background ref and elevation 
    # orientation and azimuthal angle 
    
    def __init__(self,DEM,reflectance,isWater,cell_size,alignPixels=True):
        self.DEM = DEM
        self.reflectance = reflectance
        self.isWater = isWater
        self.cell_size = cell_size
        self.alignPixels = alignPixels
        
        
        if self.DEM.shape != self.reflectance.shape:
            print('WARNING: DEM and reflectance images do not have the same shape')


        # Has to run self.set_background
        self.bg_ref = None
        self.bg_isWater = None
        self.bg_elevation = None         
        self.bg_coords = None

        # Triangulated DEM, two sets of three dimensional triangles 
        self.DEM_triangulated = None
        
        self.set_background() # test March 12, 2022
        self._triangulate_DEM()
        
        
        
    def _info(self): # print all the information 
    
        # to be continued...
    
        print("\n------- Displaying Information -------")
    
        print("\n=== DEM ")
        print(self.DEM)
        
        print("\n=== reflectance ")
        print(self.reflectance)
        
        print("\n=== isWater ")
        print(self.isWater)
        
        print("\n=== cell_size ")
        print(self.cell_size)
        
        print("\n=== bg_coords ")
        print(self.bg_coords)
        
        print("\n=== bg_elevation ")
        print(self.bg_elevation)
        
        print("\n=== bg_ref ")
        print(self.bg_ref)
        
        print("\n=== bg_isWater ")
        print(self.bg_isWater)
        
        # print("\n=== DEM_triangulated ")
        # print(self.DEM_triangulated)
        
        
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

        # default: average reflectance, average elevation, no water 
  
        
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
            self.bg_elevation = np.average(self.DEM) # default average elevation
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

        












