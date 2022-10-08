# This file is part of TMart.
#
# Copyright 2022 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.






# Replacing tmart_class

# TMart: Topography-adjusted Monte-Carlo Adjacency-effect Radiative Transfer code

from multiprocessing import cpu_count
from pathos.multiprocessing import ProcessingPool
import numpy as np
import time
import sys

# tmart dependencies 
from .tm_geometry import dirP_to_coord 
from .tm_intersect import intersect_line_DEMtri2
from .tm_water import find_R_wc, RefraIdx
from .Tmart2 import Tmart2


# Track progress in multiprocessing
def _track_job(job, update_interval=2):
    while job._number_left > 0:
        print("Tasks remaining = {0}".format(job._number_left * job._chunksize))
        time.sleep(update_interval)


# The main object in TMart
class Tmart(Tmart2):
    '''Create a Tmart object that does radiative transfer modelling. 
    
    Arguments:

    * ``Surface`` -- Surface object from the Surface module.
    * ``Atmosphere`` -- Atmosphere object from the Atmosphere module.
    * ``shadow`` -- True or False, whether to test if light is blocked by another surface at a collision. 
    * ``VROOM`` -- VROOM acceleration. 0: no acceleration. 1: all collisions are directed towards the sun. 

    Example usage::

      my_tmart = tmart.Tmart(Surface = my_surface, Atmosphere= my_atm)

    '''
    
    def __init__(self, Surface, Atmosphere, shadow = False, VROOM = 0):

        self.Surface = Surface
        self.Atmosphere = Atmosphere
        self.shadow = shadow
        self.VROOM = VROOM 
        
        self.sensor_coords = None
        self.sun_dir = None
        self.print_on = False # print switch 
        self.plot_on = False  # don't turn it on for multiprocessing 
        
        self.target_pt_direction = [180,0] 
        self.pixel = None
        self.pixel_elevation = None
        
        # Atmosphere
        self.wl = None
        self.atm_profile_wl = None # single wavelength 
        self.aerosol_SPF_wl = None 
        
        # Wind
        self.wind_speed = 10 # default 10 m/s
        self.wind_azi_avg = True # azimuthally averaged cox munk
        self.wind_dir = 0 # default azimuthal, 0 means UPWIND along x-axis
            # AKA the direction where wind comes from 
           
        # Water
        self.F_wc_wl = None # fraction of sea surface covered by whitecaps 
        self.R_wc_wl = None # whitecap reflectance at this wavelength
        
        self.water_salinity = 0
        self.water_temperature = 25      
        self.water_refraIdx_wl = None # refractive index of water at this wavelength 
        
        # In development 
        self.output_flux = False # output irradiance reflectance, direct irradiance and diffuse irradiance on the ground, under development 
        
        
             
        
    def set_geometry(self, sun_dir=[0,0], target_pt_direction=None, sensor_coords=None, pixel=None, target_coords=None):
        '''Set the sun-target-sensor geometry. Only one of ``sensor_coords``, ``pixel`` and ``target_coords`` is needed. 
        
        Arguments:
            
        * ``sun_dir`` -- Solar angle, in [Zenith, Azimuth], relative to the target.
        * ``target_pt_direction`` -- Photon's initial moving direction', AKA viewing angle, in [Zenith, Azimuth], relative to the sensor. 
            * Alternative input: 'lambertian_up'. This can only be used with ``sensor_coords``. The initial direction will be a random direction up following a Lambertian distribution for each of the photons. This is used in calculating irradiance.  
        * ``sensor_coords`` -- Where the sensor is, in [X, Y, Z], unit in meters.
        * ``pixel`` -- The target pixel to shoot photons. Parallel light rays will hit random points within the square pixel. 
        * ``target_coords`` -- Where the photon will land on the surface, only [X, Y] needed, [Z] is automatically adjusted. 


        Example usage::

          my_tmart.set_geometry(sensor_coords=[51,50,130_000], 
                                target_pt_direction=[180,0],
                                sun_dir=[0,0])
          
          # sensor_coords can be replaced by pixel or target_coords
          
          my_tmart.set_geometry(target_pt_direction=[170,0],
                                target_coords=[15,10], 
                                sun_dir=[0,0])        

          my_tmart.set_geometry(target_pt_direction=[170,0],
                                pixel=[1,1], 
                                sun_dir=[0,0])     


        '''
        
        
        n_not_none = (sum(x is not None for x in [sensor_coords,pixel,target_coords]))
        
        if n_not_none != 1: sys.exit('only one of sensor_coords,pixel,target_coords should be provided')
        
        
        
        ### Sensor coordinates, take one of the three inputs 
        
        # direct input 
        if sensor_coords is not None:
            # make sure not to hit the boundary between triangles 
            if sensor_coords[0]==sensor_coords[1]: sensor_coords[1]=sensor_coords[1]+0.0001
            self.sensor_coords = np.array(sensor_coords)            
            
            
        # pixel based, assume height 120km 
        elif pixel is not None:          
            self.pixel_elevation = self.Surface.DEM[pixel[0],pixel[1]]
            
            # distance from target to sensor 
            # it's negative because target_pt_direction is larger than 90
            dist_120000 = (120_000 - self.pixel_elevation) / np.cos(target_pt_direction[0]/180*np.pi) 
            self.sensor_coords = dirP_to_coord(dist_120000, target_pt_direction)
            self.pixel = pixel
            
            
        # target_coords, assume height 120km
        elif target_coords is not None:      
            if target_coords[0]==target_coords[1]: target_coords[1]=target_coords[1]+0.0001
            
            # target_coords is 2d, target_coords3d includes elevation
            target_coords3d = intersect_line_DEMtri2(np.array(target_coords + [120_000]), 
                                                      np.array(target_coords + [0]), 
                                                      self.Surface.DEM_triangulated)
            
            # If there is triangle intersection 
            if target_coords3d.shape[0] > 0:
                
                # closest intersection # There seems to be only 1?
                target_coords3d_chosen = target_coords3d.iloc[target_coords3d.linear_distance.idxmin()] 
                
                # print(target_coords3d)
                q_collision = target_coords3d_chosen.tolist()[0:3]  
                q_elevation = q_collision[2]
   
            # Background collision
            else:
                q_collision = np.array(target_coords + [self.Surface.bg_elevation])  
                q_elevation = q_collision[2]
            
            dist_120000 = (120_000 - q_elevation) / np.cos(target_pt_direction[0]/180*np.pi) 
            self.sensor_coords = dirP_to_coord(dist_120000, target_pt_direction) + q_collision

        else: 
            sys.exit('only one of sensor_coords,pixel,target_coords should be provided')
            
            
        # Lock photon initial direction 
        self.target_pt_direction = target_pt_direction   
        
        # Sun direction 
        self.sun_dir = sun_dir # [zenith, azimuthal]        
        
        
    def set_wind(self,wind_speed=10,wind_azi_avg=True,wind_dir=0): 
        '''Set wind speed and direction. 
        
        Arguments:

        * ``wind_speed`` -- wind speed in meters, default 10.
        * ``wind_azi_avg`` -- cox-munk slopes azimuthally averaged, default True.
        * ``wind_dir`` -- upwind direction clockwise from east, default 0. Meaningless when wind_azi_avg is True.

        Example usage::
          
          my_tmart.set_wind(wind_speed=5)
          my_tmart.set_wind(wind_speed=5, wind_azi_avg = False, wind_dir=0)
          my_tmart.set_wind(wind_speed=5, wind_azi_avg = True)
          
        '''
        
        self.wind_speed = wind_speed
        self.wind_dir = wind_dir
        self.wind_azi_avg = wind_azi_avg
        
        
    def set_water(self,water_salinity=0,water_temperature=25): # default 0/1000 and 25C
        '''Set water salinity and temperature. 
        
        Arguments:

        * ``water_salinity`` -- water salinity in parts per thousand, default 0.
        * ``water_temperature`` -- water temperature in celsius, default 25.

        Example usage::

          my_tmart.set_water(water_salinity=35, water_temperature=20)
          
        '''
        self.water_salinity = water_salinity
        self.water_temperature = water_temperature
        

    # Initiate a wavelength- or band-specific atmosphere, whitecaps properties and water refractive index 
    def _init_atm(self,band): 
        
        if self.sensor_coords is None: # Edit!!!
            print ("WARNING: geometry missing, set_geometry before you run")
        else:
            
            # Atmospheric profile and aerosol SPF
            self.atm_profile_wl, self.aerosol_SPF_wl = self.Atmosphere._wavelength(self.wl,band)
            
            # Fraction and reflectance of whitecaps 
            self.F_wc_wl, self.R_wc_wl = find_R_wc(wl=self.wl, wind_speed = self.wind_speed)
            
            # Water refractive index 
            self.water_refraIdx_wl = RefraIdx(self.water_salinity,self.water_temperature,self.wl)
            # self.water_refraIdx_wl = 1.34


    # User interface 
    def run(self, wl, band = None, n_photon=10_000, nc='auto', njobs=80, print_on=False, output_flux=False): 
        '''Run with multiple processing 
        
        Arguments:

        * ``wl`` -- wavelength in nm.
        * ``band`` -- overwrite ``wl`` with a 6S band object. We still need to specify ``wl`` because it is used in interpolating spectral SPF.
        * ``n_photon`` -- number of photons to use in MC simulation, default 10,000.
        * ``nc`` -- number of CPU cores to use in multiprocessing, default automatic. 
        * ``njobs`` -- dividing the jobs into n portions in multiprocessing, default 80. 
        
        Return:

        * Movement information of photons.
        
        Example usage::

          n_photon = 100_00
          nc = 10
          njobs = 100
          results = my_tmart.run(wl=wl, n_photon=n_photon,nc= nc,njobs= njobs)
          
        '''
        
        self.wl = wl 
        self.print_on = print_on
        self.plot_on = False # don't even try it 
        self.output_flux = output_flux
        self._init_atm(band)
        
        
        if nc=='auto':
            nc = cpu_count()
        else:
            nc = nc
            
        print("\n========= Initiating T-Mart =========")
        print(f"Number of photons: {n_photon}")
        print(f'Using {nc} core(s)')
        
        n = n_photon
        
        
        part_count = [n/njobs for i in range(njobs)]
        
        part_count = np.array_split(range(n_photon), njobs)
        
        
        print(f"Number of job(s): {njobs}")
        print('Wavelength: ' + str(self.wl))
        print('target_pt_direction: ' + str(self.target_pt_direction))
        print('sun_dir: ' + str(self.sun_dir))
        print("=====================================")
        

        pool = ProcessingPool(processes=nc)
        time.sleep(0.5)
        
        # manual print 
        results_temp = pool.amap(self._run,part_count) # Async
        
        if njobs>1:
            _track_job(results_temp)
        
        results = results_temp.get()
        results = np.vstack(results)

        return results 

        
    # Distribute runs to processors     
    def _run(self,part_count):
    
        pts_stat = np.empty([0,13])
        
        for i in part_count:
            
            if self.print_on:
                print("\n---------- Running Photon " + str(i) + " ----------")
            
            pt_stat = self._run_single_photon(i)
            pts_stat = np.vstack([pts_stat, pt_stat])
      
        return pts_stat
    
    
    def run_plot(self, wl, band = None, plot_on=True, plot_range=[0,100_000,0,100_000,0,100_000]): 
        '''Run a single photon and plot, print the details of photon movements. 
        To observe the photon movements, mostly for debugging purposes. 
        
        Arguments:

        * ``wl`` -- wavelength in nm.
        * ``band`` -- overwrite ``wl`` with a 6S band object. We still need to specify ``wl`` because it is used in interpolating spectral SPF.
        * ``plot_on`` -- Boolean, if plot the movements. 
        * ``plot_range`` -- List, [xmin, xmax, ymin, ymax, zmin, zmax]
        
        Return:

        * Movement information of a photon.

        Example usage::

          results = my_tmart.run_plot(wl=wl, plot_on=True, plot_range=[0,100_000,0,100_000,0,100_000])
          
        '''
        
        
        print("\n====== Running and Plotting T-Mart Single Photon ======")
        self.wl = wl 
        self.print_on = True    # Always print the details of photon movements 
        self.plot_on = plot_on  # Default plot, may turn off 
        self.plot_range = plot_range
        self._init_atm(band)
        
        return self._run_single_photon(0)
    
    











