# This file is part of TMart.
#
# Copyright 2022 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.




# Atmosphere object  

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

import Py6S
from Py6S.Params.atmosprofile import AtmosProfile
from Py6S.Params.aeroprofile import AeroProfile

from .Aerosol import find_aerosolSPF

import os.path


class Atmosphere(): 
    '''Create an Atmosphere object that is wavelength independent. I.e., a wavelength-dependent atmosphere is generated every time TMart is run.
    
    Arguments:

    * ``atm_profile`` -- AtmosProfile object from Py6S. 
    * ``aot550`` -- AOT at 550nm.
    * ``aerosol_type`` -- 'BiomassBurning', 'Continental', 'Desert', 'Maritime', 'Stratospheric' or 'Urban', as provided by 6S. 
    * ``wl`` -- central wavelength in nm.
    * ``n_layers`` -- Number of atmosphere layers to use. Default 20. 
    * ``AEROSOL_SCALE_HEIGHT`` -- Aerosol scale height in km. Default 2km. 
    * ``no_absorption`` -- Boolean, if yes -> remove all absorption. 
    * ``specify_ot_rayleigh`` -- Specify rayleigh optical thickness, only for testing.
    * ``specify_abs`` -- Specifiy absorption optical thickness, only for testing. 
    

    Example usage:: 
    
      from Py6S.Params.atmosprofile import AtmosProfile

      # Atmophere profile comes from 6S
      atm_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeSummer) 
      aerosol_type = 'Maritime' 
      aot550 = 0.1
      n_layers = 20
      aerosol_scale_height = 2 # Unless you have a reason, don't change this
        
      # Synthesize an atmosphere object    
      my_atm = tmart.Atmosphere(atm_profile, aot550, aerosol_type, n_layers, aerosol_scale_height)
    '''

    def __init__(self,atm_profile, aot550 = 0, 
                 aerosol_type='Maritime' , 
                 n_layers=20, AEROSOL_SCALE_HEIGHT=2, no_absorption = False, specify_ot_rayleigh = -1, specify_abs = -1):
        
        self.atm_profile = atm_profile
        self.aot550 = aot550
        
        
        self.aerosol_type = aerosol_type
        
        # Default 20 layers 
        self.n_layers = n_layers
        
        # Default aerosol scale height is 2km
        self.aerosol_scale_height = AEROSOL_SCALE_HEIGHT
     
        # atm is always 100km thick in total, following 6S 
        self.atm_height = 100
        self.layer_height = self.atm_height/self.n_layers

        # special scenarios 
        self.no_absorption = no_absorption
        self.specify_ot_rayleigh = specify_ot_rayleigh
        self.specify_abs = specify_abs
        
    
    ### The main one dealing with wavelength and band 
    
    # return everything at a single wavelength 
    # supports Aerosol LUT AND 6S wavelengths 
    # Return a table of OTs at different heights  +  aerosol_SPF
    def _wavelength(self, wl, band=None): 
        
        self.wl = wl
        
        # Find bottom, mean and top altitudes of layers 
        self.layers_alts_bottom = np.linspace(0, self.atm_height - self.layer_height, self.n_layers)  
        self.layers_alts_mean = np.linspace(self.layer_height/2, self.atm_height - self.layer_height/2, self.n_layers)  
        self.layers_alts_top = np.linspace(self.layer_height, self.atm_height , self.n_layers)  
        
        # Calculate tau in each layer and extract the molecular profile at one wavelength 
        layers_ot_molecule, layers_ot_rayleigh = self._atm_profile_wl(band)
        
        # ot_mie is aerosol scattering, ot_aerosol is aerosol absorption 
        layers_ot_mie, layers_ot_aerosol = self._aerosol_wl(band)

        # Retrieve aerosol SPF
        aerosol_SPF = find_aerosolSPF(self.aerosol_type,wl)
    
        layers_ot_total_abs = layers_ot_molecule + layers_ot_aerosol
        
        # A data frame of optical thickness 
        atm_OT = pd.DataFrame({'Alt_bottom': self.layers_alts_bottom,
                            'Alt_top': self.layers_alts_top,  
                            'ot_abs': layers_ot_total_abs, 
                            'ot_rayleigh': layers_ot_rayleigh,
                            'ot_mie': layers_ot_mie
                            })
        
        if self.specify_ot_rayleigh == -1:
            pass
        else:
            atm_OT['ot_rayleigh'] = self.specify_ot_rayleigh       
        
        
        atm_OT['ot_scatt'] = atm_OT.ot_rayleigh + atm_OT.ot_mie
        atm_OT['l_height'] = atm_OT.Alt_top - atm_OT.Alt_bottom
        atm_OT['percentage'] = 0 # used to calculate travelling distance 
        
        
        if self.no_absorption:
            atm_OT['ot_abs'] = 0
            
        if self.specify_abs == -1:
            pass
        else:
            atm_OT['ot_abs'] = self.specify_abs               
        
        return atm_OT, aerosol_SPF
    
    
    
    # Extract the molecular profile at one wavelength 
    # Return two lists: ot_molecule and ot_rayleigh 
    def _atm_profile_wl(self,band): 

        layers_ot_molecule = []
        layers_ot_rayleigh = []
        
        layers_alts_bottom = self.layers_alts_bottom
        layers_alts_top = self.layers_alts_top
           
        
        for i_layer in range(len(layers_alts_bottom)):
        
            alt_bottom = 0
            alt_top = layers_alts_top[i_layer]
            
            # Initiate an SixS object 
            s = Py6S.SixS()
            
            # Wavelength and band 
            s.wavelength = Py6S.Wavelength(self.wl/1000)   
            if band is not None:
                s.wavelength = band
            
            # Other parameters 
            s.atmos_profile = self.atm_profile        
            s.geometry = Py6S.Geometry.User()
            s.geometry.solar_z, s.geometry.solar_a = 0, 0
            s.geometry.view_z, s.geometry.view_a = 0, 0
            s.aero_profile = AeroProfile.PredefinedType(AeroProfile.NoAerosols)
            
            # Altitudes 
            s.altitudes.set_target_custom_altitude(alt_bottom)
            
            if alt_top < 100:
                s.altitudes.set_sensor_custom_altitude(alt_top)
            else:
                s.altitudes.set_sensor_satellite_level()    
                s.altitudes.set_sensor_custom_altitude(99.99999)
        
            s.run()
            
            # print(s.outputs.fulltext)
            # print(s.outputs.optical_depth_total)
            # print(s.outputs.optical_depth_total.rayleigh)
    
            tao_rayleigh = s.outputs.optical_depth_plane.rayleigh
    
            # Molecular absorption 
            T_molecule = s.outputs.transmittance_global_gas.upward 
            tao_molecule = -np.log(T_molecule)
            
            # we keep last tao_molecule and tao_rayleigh for full atmosphere 
            layers_ot_molecule.append(tao_molecule)
            layers_ot_rayleigh.append(tao_rayleigh)
                
        
        # Convert 0-10, 0-20... to 0-10, 10-20...
        layers_ot_molecule_new = [layers_ot_molecule[0]]
        layers_ot_rayleigh_new = [layers_ot_rayleigh[0]]    
        
        for i_layer in range(1, len(layers_ot_molecule)):
            
            i_molecule = layers_ot_molecule[i_layer] - layers_ot_molecule[i_layer-1] 

            if i_molecule < 0:
                i_molecule = 0
            
            layers_ot_molecule_new.append(i_molecule)
            
            i_rayleigh = layers_ot_rayleigh[i_layer] - layers_ot_rayleigh[i_layer-1] 

            if i_rayleigh < 0:
                i_rayleigh = 0
            
            layers_ot_rayleigh_new.append(i_rayleigh)            
            
        # percentage in each layer, redistribute total  
        if sum(layers_ot_molecule_new)> 0: # no action in case of full transparency 
            relative_mol = layers_ot_molecule_new/np.sum(layers_ot_molecule_new)
            layers_ot_molecule_new = tao_molecule * relative_mol  
        if sum(layers_ot_rayleigh_new)> 0: 
            relative_ray = layers_ot_rayleigh_new/np.sum(layers_ot_rayleigh_new)    
            layers_ot_rayleigh_new = tao_rayleigh * relative_ray


        return layers_ot_molecule_new, layers_ot_rayleigh_new   
    
    


    # Extract aerosol profile at one wavelength 
    # Find the spectral dependence of AOT and then apply it to AOT550
    def _aerosol_wl(self, band):
        
        # Initiate an SixS object
        s = Py6S.SixS()
        
        # Band and wavelength 
        s.wavelength = Py6S.Wavelength(self.wl/1000)   
        if band is not None:
            s.wavelength = band
        
        # Other parameters 
        s.atmos_profile = self.atm_profile        
        s.geometry = Py6S.Geometry.User()
        s.geometry.solar_z, s.geometry.solar_a = 0, 0
        s.geometry.view_z, s.geometry.view_a = 0, 0
        s.altitudes.set_target_custom_altitude(0)
        s.altitudes.set_sensor_satellite_level()
        
        aerosol_dict = {    "NoAerosols": 0,
                        "Continental" : 1,
                        "Maritime" : 2,
                        "Urban" : 3,
                        "Desert" : 5,
                        "BiomassBurning" : 6,
                        "Stratospheric" : 7}
        
        
        s.aero_profile = AeroProfile.PredefinedType(aerosol_dict[self.aerosol_type])
        s.aot550 = 1
        s.run()
        
    
        aerosol_EXT_r550 = s.outputs.optical_depth_total.aerosol # ratio, relative to 550
        
        aerosol_EXT = self.aot550 * aerosol_EXT_r550 # total extinction 
        aerosol_SSA = s.outputs.single_scattering_albedo.aerosol # single scattering albedo 
        
        ot_mie = aerosol_EXT * aerosol_SSA # scattering by aerosols 
        ot_aerosol = aerosol_EXT - ot_mie  # absorption by aerosols 
        
        # Concentration at mean altitudes, make integrated area???
        conc_relative = np.exp(-self.layers_alts_mean/self.aerosol_scale_height)
        conc_normalized = conc_relative / np.sum(conc_relative)        

        # Divide optical thickness into layers, weighted by concentration at mean heights 
        layers_ot_mie = ot_mie * conc_normalized
        layers_ot_aerosol = ot_aerosol * conc_normalized # aerosol absorption
  
        return layers_ot_mie, layers_ot_aerosol # scattering and absorption by aerosol, respectively 











