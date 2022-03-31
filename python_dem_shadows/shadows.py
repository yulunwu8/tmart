import numpy as np
import sys


def project_shadows(dem, sun_vector, dx, dy=None, print_progress=False):
    """Cast shadows on the DEM from a given sun position."""

    if dy is None:
        dy = dx

    inverse_sun_vector = _invert_sun_vector(sun_vector)
    normal_sun_vector = _normalize_sun_vector(sun_vector)

    rows, cols = dem.shape # n_row, n_cols
    z = dem.T

    # Determine sun direction.
    if sun_vector[0] < 0:
        # The sun shines from the West.
        start_col = 1
    else:
        # THe sun shines from the East.
        start_col = cols - 1

    if sun_vector[1] < 0:
        # The sun shines from the North.
        start_row = 1
    else:
        # The sun shines from the South.
        start_row = rows - 1

    in_sun = np.ones_like(z)
    
    
    # if the lowest point of the DEM is 0, we speed up the computation by skipping water's shadows 
    if np.min(dem) == 0: 
        initiate = True 
    else:
        initiate = False
    
    
    # Project North-South
    row = start_row
    for col in range(cols):
        if print_progress: print('Project North-South: {}/{}'.format(col, cols))
        _cast_shadow(row, col, rows, cols, dx, in_sun, inverse_sun_vector,
                     normal_sun_vector, z, initiate)

    # Project West-East
    col = start_col
    for row in range(rows):
        if print_progress: print('Project West-East: {}/{}'.format(row, rows))
        _cast_shadow(row, col, rows, cols, dy, in_sun, inverse_sun_vector,
                     normal_sun_vector, z, initiate)
    return in_sun.T


def _normalize_sun_vector(sun_vector):
    normal_sun_vector = np.zeros(3)
    normal_sun_vector[2] = np.sqrt(sun_vector[0] ** 2 + sun_vector[1] ** 2)
    normal_sun_vector[0] = -sun_vector[0] * sun_vector[2] / normal_sun_vector[2]
    normal_sun_vector[1] = -sun_vector[1] * sun_vector[2] / normal_sun_vector[2]
    return normal_sun_vector


def _invert_sun_vector(sun_vector):
    return -sun_vector / max(abs(sun_vector[:2]))


def _cast_shadow(row, col, rows, cols, dl, in_sun, inverse_sun_vector,
                 normal_sun_vector, z, initiate):
    n = 0
    z_previous = -sys.float_info.max
    """
    Explanation Hans
    ------------
    Inverse sun vector allows you to calculate the projection on projection plane back to the original matrix.

    The dx and dy represent how far the shadow extends from the current cell. 
    Normal sun vector is the plane perpendicular to the direction of the sun.
    """
    
    
    # 'initiate' added by Yulun Wu: to skip the elevations 0s in the begining of the line, as
    # often encountered in coastal areas. 
    _initiate = True
    

    while True:
        # Calculate projection offset
        dx = inverse_sun_vector[0] * n
        dy = inverse_sun_vector[1] * n
        col_dx = int(round(col + dx))
        row_dy = int(round(row + dy))
        if (col_dx < 0) or (col_dx >= cols) or (row_dy < 0) or (row_dy >= rows):
            break
        

        vector_to_origin = np.zeros(3)
        vector_to_origin[0] = dx * dl
        vector_to_origin[1] = dy * dl
        vector_to_origin[2] = z[col_dx, row_dy]

        if not _initiate:
            if vector_to_origin[2]>0:
                _initiate = True
        else:
            z_projection = np.dot(vector_to_origin, normal_sun_vector)
            
            # YW understanding: z_projection is the projected height of the sun 

            if z_projection < z_previous:
                in_sun[col_dx, row_dy] = 0
            else:
                z_previous = z_projection
                
        n += 1








