# This file is part of TMart.
#
# Copyright 2024 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# Source: https://github.com/acolite/acolite/blob/main/acolite/shared/fillnan.py

def fillnan(data):
    from scipy.ndimage import distance_transform_edt
    import numpy as np

    ## fill nans with closest value
    ind = distance_transform_edt(np.isnan(data), return_distances=False, return_indices=True)
    return(data[tuple(ind)])
