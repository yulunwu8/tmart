# Atmosphere object  

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

import Py6S
from Py6S.Params.atmosprofile import AtmosProfile
from Py6S.Params.aeroprofile import AeroProfile

from .Aerosol import find_aerosolSPF

import os.path

# Include layered atmosphere, molecules and particles, at all wavelengths  

class Atmosphere(): # wavelength
    '''Create an Atmosphere object. 
    
    Arguments:

    * ``atm_profile`` -- AtmosProfile object from Py6S.
    * ``aot550`` -- AOT at 550nm.
    * ``aerosol_type`` -- 'BiomassBurning', 'Continental', 'Desert', 'Maritime', 'Stratospheric' or 'Urban', as provided by 6S. 
    * ``wl`` -- central wavelength in nm.
    * ``n_layers`` -- Number of atmosphere layers to use
    * ``aerosol_scale_height`` -- Aerosol scale height in km. Default 2km. 
    * ``no_absorption`` -- Boolean, if yes -> remove all absorption. 
    * ``specify_ot_rayleigh`` -- Specify rayleigh optical thickness.
    * ``specify_abs`` -- Specifiy absorption optical thickness.
    

    Example usage:: ### EDIT!!!
    
      from Py6S.Params.atmosprofile import AtmosProfile

      atm_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeSummer) 
      aerosol_SPF = 'tmart/ancillary/aerosol_maritime_SPF.csv' 
      aot550 = 0.1
      my_atm = tmart.Atmosphere(atm_profile, aot550, aerosol_SPF)

    '''
    
    # os.path.join(os.path.dirname(__file__), 'ancillary/aerosol_maritime_SPF.csv')

    def __init__(self,atm_profile, aot550 = 0, 
                 aerosol_type='Maritime' , wl = None, 
                 n_layers=None, aerosol_scale_height=None, no_absorption = False, specify_ot_rayleigh = -1, specify_abs = -1):
        
        self.atm_profile = atm_profile
        self.aot550 = aot550
        
        # pick a name and these are loaded automatically 
        
        self.aerosol_type = aerosol_type
        
        
        # SPF is from HERE!!! 
        
        self.aerosol_SPF = find_aerosolSPF(aerosol_type,wl)
        
        # aerosol_SPF = os.path.join(os.path.dirname(__file__), 'ancillary/aerosolSPF_maritime.csv')
        # self.aerosol_SPF = pd.read_csv(aerosol_SPF) 

        
        # Default 20 layers 
        if n_layers is None:
            self.n_layers = 20
        else:
            self.n_layers = n_layers
        
        # Default aerosol scale height is 2km
        if aerosol_scale_height is None:
            self.aerosol_scale_height = 2
        else:
            self.aerosol_scale_height = aerosol_scale_height        
        
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
    # A table of OTs at different heights   +    aerosol_SPF
    def _wavelength(self, wl, band=None): 
        
        self.wl = wl
        
        # Find bottom, mean and top altitudes of layers 
        self.layers_alts_bottom = np.linspace(0, self.atm_height - self.layer_height, self.n_layers)  
        self.layers_alts_mean = np.linspace(self.layer_height/2, self.atm_height - self.layer_height/2, self.n_layers)  
        self.layers_alts_top = np.linspace(self.layer_height, self.atm_height , self.n_layers)  
        
        # print (self.layers_alts_bottom)
        # print (self.layers_alts_mean)
        # print (self.layers_alts_top)
        # print (self.layers_alts_top[-1]) # print last one 
        # print (type(self.layers_alts_top))
        
        
        
        # test = self._aerosol_wl()
        # print(test[0])
        # print(test[1])
        
        # print(sum(test[0]))
        # print(sum(test[1]))
        # print(sum(test[0])+ sum(test[1]))
        
        # Calculate tau in each layer and extract the molecular profile at one wavelength 
        layers_ot_molecule, layers_ot_rayleigh = self._atm_profile_wl(band)
        
        # print('\n======= Summary ======')
        # print(layers_ot_molecule)
        # print(layers_ot_rayleigh)
        # print(sum(layers_ot_molecule))
        # print(sum(layers_ot_rayleigh))
        
        
        # ot_mie is aerosol scattering, ot_aerosol is aerosol absorption 
        layers_ot_mie, layers_ot_aerosol, aerosol_SPF = self._aerosol_wl(band)
        # print(layers_ot_aerosol)
    
        layers_ot_total_abs = layers_ot_molecule + layers_ot_aerosol
        
        # print(self.layers_alts_bottom)
        # print(self.layers_alts_top)
        # print(layers_ot_total_abs)
        # print(layers_ot_rayleigh)
        # print(layers_ot_mie)
        
        atm_OT = pd.DataFrame({'Alt_bottom': self.layers_alts_bottom,
                            'Alt_top': self.layers_alts_top,  
                            'ot_abs': layers_ot_total_abs, 
                            'ot_rayleigh': layers_ot_rayleigh,
                            'ot_mie': layers_ot_mie
                            })
        
        # # Customize atm_OT
        # atm_OT = pd.DataFrame({'Alt_bottom': [0],
        #                     'Alt_top': [100],  
        #                     'ot_abs': [0], 
        #                     'ot_rayleigh': [0.1],
        #                     'ot_mie': [0.1]
        #                     })        
        
        # atm_OT.ot_abs = 0.0
        # atm_OT.ot_rayleigh = 0.143263
        
        # atm_OT = pd.read_csv('atm_lRt400_noAbs.csv')
        
        
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
            
            # print('alt_bottom: ' + str(alt_bottom))
            # print('alt_top: ' + str(alt_top))
            
            s = Py6S.SixS()
            s.wavelength = Py6S.Wavelength(self.wl/1000)   
            
            
            ### HERE
            
            if band is not None:
                s.wavelength = band
            
            # s.wavelength = Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_08)
            
            
            
            s.atmos_profile = self.atm_profile        
        
            s.geometry = Py6S.Geometry.User()
            s.geometry.solar_z, s.geometry.solar_a = 0, 0
            s.geometry.view_z, s.geometry.view_a = 0, 0
            s.aero_profile = AeroProfile.PredefinedType(AeroProfile.NoAerosols)
            
            # altitudes 
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
            # print('tao_rayleigh: ' + str(tao_rayleigh))
    
    
            # Molecular absorption 
            T_molecule = s.outputs.transmittance_global_gas.upward 
            # print('T_molecule: ' + str(T_molecule))
            
            tao_molecule = -np.log(T_molecule)
            
            # print('tao_molecule: ' + str(tao_molecule))
            
            # we keep last tao_molecule and tao_rayleigh for full atmosphere 
            
            layers_ot_molecule.append(tao_molecule)
            layers_ot_rayleigh.append(tao_rayleigh)
            
        # print('layers_ot_molecule: ' + str(layers_ot_molecule))
        # print('layers_ot_rayleigh: ' + str(layers_ot_rayleigh))
        
        
        # Convert 0-10, 0-20... to 0-10, 10-20...
        layers_ot_molecule_new = [layers_ot_molecule[0]]
        layers_ot_rayleigh_new = [layers_ot_rayleigh[0]]    
        
        for i_layer in range(1, len(layers_ot_molecule)):
            
            # print('==========')
            # print(i_layer)
            
            # print(layers_ot_molecule[i_layer])
            # print(layers_ot_molecule[i_layer-1]) 
            
            
            i_molecule = layers_ot_molecule[i_layer] - layers_ot_molecule[i_layer-1] 
            # print(i_molecule)
            if i_molecule < 0:
                i_molecule = 0
            
            layers_ot_molecule_new.append(i_molecule)
            
            
            i_rayleigh = layers_ot_rayleigh[i_layer] - layers_ot_rayleigh[i_layer-1] 
            # print(i_rayleigh)
            if i_rayleigh < 0:
                i_rayleigh = 0
            
            layers_ot_rayleigh_new.append(i_rayleigh)            
            
        
        
        # percentage in each layer, redistribute total  
        
        if sum(layers_ot_molecule_new)> 0: # no action in case of full transparency 
            relative_mol = layers_ot_molecule_new/np.sum(layers_ot_molecule_new)
            
            # tao_molecule = 0.3
            
            layers_ot_molecule_new = tao_molecule * relative_mol  
        
        
        if sum(layers_ot_rayleigh_new)> 0: 
            relative_ray = layers_ot_rayleigh_new/np.sum(layers_ot_rayleigh_new)
            
            # test
            # tao_rayleigh = 0.3601303
            
            layers_ot_rayleigh_new = tao_rayleigh * relative_ray
            
        # print('layers_ot_molecule_new: ' + str(layers_ot_molecule_new))
        # print(sum(layers_ot_molecule_new))
        # print('layers_ot_rayleigh_new: ' + str(layers_ot_rayleigh_new))
        # print(sum(layers_ot_rayleigh_new))

        return layers_ot_molecule_new, layers_ot_rayleigh_new   
    
    


    # Extract aerosol profile at one wavelength 
    def _aerosol_wl(self, band):
        # Adapt to aerosol model???
        
        
        s = Py6S.SixS()
        s.wavelength = Py6S.Wavelength(self.wl/1000)   
        
        
        # HERE
        
        # s.wavelength = Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_08)
        
        if band is not None:
            s.wavelength = band
        
        
        
        s.atmos_profile = self.atm_profile        
    
        s.geometry = Py6S.Geometry.User()
        s.geometry.solar_z, s.geometry.solar_a = 0, 0
        s.geometry.view_z, s.geometry.view_a = 0, 0
        # s.aero_profile = AeroProfile.PredefinedType(AeroProfile.NoAerosols)
        
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
        
    
        aerosol_EXT_r550 = s.outputs.optical_depth_total.aerosol # relative to 550
        
        aerosol_EXT = self.aot550 * aerosol_EXT_r550 # total extinction 
        aerosol_SSA = s.outputs.single_scattering_albedo.aerosol
        
        # print('aerosol_EXT_r550: ' + str(aerosol_EXT_r550))
        # print('aerosol_SSA: ' + str(aerosol_SSA))
    
    
        ### Old Method 
        
        # # Interpolating extinction coefficient 
        # ext_wl = self.aerosol_EXT.Wavelength.to_numpy()
        # ext_value = self.aerosol_EXT.Value.to_numpy()
        # ext_f = interp1d(ext_wl, ext_value, kind='cubic')
        
        # # Interpolating single-scattering albedo 
        # ssa_wl = self.aerosol_SSA.Wavelength.to_numpy()
        # ssa_value = self.aerosol_SSA.Value.to_numpy()
        # ssa_f = interp1d(ssa_wl, ssa_value, kind='cubic')        
        
        # # Calculations 
        # aerosol_EXT_r550 = ext_f(self.wl).item()     
        # aerosol_EXT = self.aot550 * aerosol_EXT_r550
        # aerosol_SSA = ssa_f(self.wl).item()
        

        
        
        ot_mie = aerosol_EXT * aerosol_SSA # scattering by aerosols 
        ot_aerosol = aerosol_EXT - ot_mie  # absorption by aerosols 
        
        # Concentration at mean altitudes, make integrated area???
        conc_relative = np.exp(-self.layers_alts_mean/self.aerosol_scale_height)
        conc_normalized = conc_relative / np.sum(conc_relative)        
        
        # print(conc_relative)
        # print(conc_normalized)
        # print(np.sum(conc_normalized))

        # Divide optical thickness into layers, weighted by concentration at mean heights 
        layers_ot_mie = ot_mie * conc_normalized
        layers_ot_aerosol = ot_aerosol * conc_normalized # aerosol absorption

        # Normalized scattering phase function
        aerosol_SPF = self.aerosol_SPF
        
        
        # print('layers_ot_mie: ' + str(layers_ot_mie))
        # print('layers_ot_aerosol: ' + str(layers_ot_aerosol))
        # print(sum(layers_ot_aerosol))
  
    
  
        return layers_ot_mie, layers_ot_aerosol, aerosol_SPF



if __name__=='__main__':


    atm_profile = AtmosProfile.PredefinedType(AtmosProfile.MidlatitudeSummer) # same as 6S
    aot550 = 0.0
    # aerosol_SPF = 'aerosol_maritime_SPF.csv'
    # aerosol_EXT = 'aerosol_maritime_EXT.csv'
    # aerosol_SSA = 'aerosol_maritime_SSA.csv'
    

    # all wavelengths 
    my_atm = Atmosphere(atm_profile, aot550)
    
    
    # my_atm.wavelength(wl=550,n_layers=10)
        
    '''
    
    my_atm.wavelength(wl=500)
    

    
    '''
    
    band = Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_08)
    


    my_atm_OT, my_SPF  = my_atm._wavelength(wl=842, band = band)
    
    print(sum(my_atm_OT['ot_abs']))
    print(sum(my_atm_OT['ot_rayleigh']))
    print(sum(my_atm_OT['ot_scatt']))
    
    
    # my_atm_OT_band, my_SPF_band  = my_atm._wavelength(wl=842) 
    
    # sum(my_atm_OT_band['ot_abs'])
    # sum(my_atm_OT_band['ot_rayleigh'])
    # sum(my_atm_OT_band['ot_scatt'])    
    
    
    # my_atm_OT2, my_SPF2  = my_atm.wavelength(wl=950)
    
    # my_atm_OT2, my_SPF2  = my_atm.wavelength(wl=550,aerosol_scale_height=10, n_layers=20)
    
    # my_atm_OT.to_csv('test_atm_profile2.csv', index=False)
    
    # print(sum(my_atm_OT.ot_abs))
    

    
    









