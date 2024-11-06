# This file is part of T-Mart.
#
# Copyright 2024 Yulun Wu.
#
# T-Mart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# Compute total gas transmittance 

from Py6S import *
import numpy as np

def normal_distribution(x, mu, sigma):
    return 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-((x - mu)**2) / (2 * sigma**2))

def calculate_heights(center_wavelength, FWHM, range_width=100, interval=2.5):
    sigma = FWHM / (2 * np.sqrt(2 * np.log(2)))
    start_wavelength = center_wavelength - range_width / 2
    end_wavelength = center_wavelength + range_width / 2
    wavelengths = np.arange(start_wavelength, end_wavelength + interval, interval)
    heights = normal_distribution(wavelengths, center_wavelength, sigma)
    
    # normalize heights 
    heights = heights / np.max(heights)
    heights = np.round(heights, 6)

    return wavelengths[0], wavelengths[-1], heights

def compute_gas_transmittance(metadata, anci, wavelength, FWHM):

    s = SixS()
    s.geometry = Geometry.User()
    s.geometry.solar_z = metadata['sza']
    s.geometry.solar_a = 0
    s.geometry.view_z = 0 
    s.geometry.view_a = 0
    s.atmos_profile = AtmosProfile.UserWaterAndOzone(anci['water_vapour']/10, anci['ozone']/1000)
    s.altitudes.set_sensor_satellite_level()
    wl_start, wl_end, heights = calculate_heights(wavelength,FWHM)
    s.wavelength = Wavelength(wl_start/1000, wl_end/1000, heights.tolist())
    
    s.run()

    return s.outputs.transmittance_global_gas.total
