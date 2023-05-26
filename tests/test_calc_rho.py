# This file is part of TMart.
#
# Copyright 2022 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# %matplotlib qt

# Go up by 2 directory and import 

import sys
import os.path as path
two_up =  path.abspath(path.join(__file__ ,"../.."))
sys.path.append(two_up)


# Actual imports 
import tmart
import numpy as np
from Py6S.Params.atmosprofile import AtmosProfile
import Py6S


wl = [400,800,100]

wl = 800


rho = tmart.surface_rho.calculate(wl=wl, viewing_zenith=30, solar_zenith=30.88, 
                                  relative_azimuth=90, aot550=0.05, wind_speed = 8, n_photon=10_000)


print(rho)






