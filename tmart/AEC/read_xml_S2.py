# This file is part of TMart.
#
# Copyright 2024 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# Read S2 metadata in XML

def read_xml_S2(file):
    
    from xml.dom import minidom
    
    # Parse an xml file by name
    xml = minidom.parse(file)
    
    # General
    General_Info = xml.getElementsByTagName('n1:General_Info')
    time = General_Info[0].getElementsByTagName('SENSING_TIME')[0].firstChild.nodeValue

    # Geometric 
    Geometric_Info = xml.getElementsByTagName('n1:Geometric_Info')

    for ta in Geometric_Info[0].getElementsByTagName('Tile_Angles'):
        for sub in ta.getElementsByTagName('Mean_Sun_Angle'):
            mean_sza = float(sub.getElementsByTagName('ZENITH_ANGLE')[0].firstChild.nodeValue)
            mean_saa = float(sub.getElementsByTagName('AZIMUTH_ANGLE')[0].firstChild.nodeValue)
        vzas = []
        vaas = []
        for sub in ta.getElementsByTagName('Mean_Viewing_Incidence_Angle_List'):
            for band in ta.getElementsByTagName('Mean_Viewing_Incidence_Angle'):
                vza = float(band.getElementsByTagName('ZENITH_ANGLE')[0].firstChild.nodeValue)
                vaa = float(band.getElementsByTagName('AZIMUTH_ANGLE')[0].firstChild.nodeValue)
                vzas.append(vza)
                vaas.append(vaa)
                
    mean_vza = sum(vzas) / len(vzas)
    mean_vaa = sum(vaas) / len(vaas)
    tm_pt_dir = [180-mean_vza, (mean_vaa+90)%360] # mean_vaa = 0 => we are looking north, photon move south => 90 in T-Mart
    tm_sun_dir=[mean_sza, (mean_saa+270)%360] # mean_saa = 0 => sun is in the north => 270 in T-Mart
    
    return {'time': time,
            'sza': mean_sza,
            'saa': mean_saa,
            'vza': mean_vza,
            'vaa': mean_vaa,
            'tm_pt_dir': tm_pt_dir,
            'tm_sun_dir': tm_sun_dir}

