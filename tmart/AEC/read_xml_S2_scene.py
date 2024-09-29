# This file is part of TMart.
#
# Copyright 2024 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# Read S2 metadata in XML

def read_xml_S2_scene(file):
    
    from xml.dom import minidom
    import sys
    
    # Parse an xml file by name
    xml = minidom.parse(file)
    bands = ['B01','B02','B03','B04','B05','B06','B07','B08','B8A','B09','B10','B11','B12']
    metadata = {}
    tdom = xml.getElementsByTagName('RADIO_ADD_OFFSET')
    
    # Older imagery 
    if len(tdom) == 0:
        for band_name in bands:
            metadata[str(band_name) + '_mult'] = 1/10_000
            metadata[str(band_name) + '_add'] = 0
    
    else: 
        for i in range(len(tdom)):
            sub_tdom = tdom[i]
            band_id = int(sub_tdom.getAttribute('band_id'))
            if band_id == i:
                band_name = bands[i]
                metadata[str(band_name) + '_mult'] = 1/10_000
                metadata[str(band_name) + '_add'] = float(sub_tdom.firstChild.nodeValue) / 10_000
            else: 
                sys.exit('Warning: failed to read {}'.format(file))
    return metadata   
            