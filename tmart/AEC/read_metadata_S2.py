# This file is part of TMart.
#
# Copyright 2024 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# Read S2 metadata
# With contribution from Shun Bi

def read_metadata_S2(file,config,sensor):
    import tmart
    import os, sys
    import mgrs
    import Py6S, math
    
    # search metadata xml file
    for root, dirs, files in os.walk(file):
        if "MTD_MSIL1C.xml" in files:
            xml_file = os.path.join(root, "MTD_MSIL1C.xml")
    
    file = os.path.dirname(xml_file)
    files = os.listdir(file)
    metadata = {}
    metadata['file'] = file
    
    # resolution for masks 
    mask_res = [10,20,60] # mask resolution 
    mask_res.append(int( 10 *  int(config['reshape_factor_S2'])))
    metadata['mask_res'] = mask_res
    
    # high TOA bands for masks 
    highTOA_band_names = ['B01','B02','B03','B04','B05','B06','B07','B08','B8A','B11','B12']
    highTOA_band_names.remove(config['S2_SWIR_band'])
    metadata['highTOA_band_names'] = highTOA_band_names
    
    # masks 
    metadata['cirrus_mask'] = config['S2_cirrus_band']
    metadata['SWIR_mask'] = config['S2_SWIR_band']
    
    # bands to be AECed
    metadata['AEC_bands_name'] = ['B01','B02','B03','B04','B05','B06','B07','B08','B8A','B11','B12']
    
    if sensor == 'S2A': 
        metadata['AEC_bands_6S'] = [Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_01),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_02),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_03),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_04),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_05),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_06),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_07),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_08),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_8A),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_11),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2A_MSI_12)]
        metadata['AEC_bands_wl'] = [442.7, 492.7, 559.8, 664.6, 704.1, 740.5, 782.8, 832.8, 864.7, 1613.7, 2202.4]
    elif sensor == 'S2B':
        metadata['AEC_bands_6S'] = [Py6S.Wavelength(Py6S.PredefinedWavelengths.S2B_MSI_01),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2B_MSI_02),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2B_MSI_03),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2B_MSI_04),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2B_MSI_05),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2B_MSI_06),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2B_MSI_07),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2B_MSI_08),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2B_MSI_8A),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2B_MSI_11),
                                    Py6S.Wavelength(Py6S.PredefinedWavelengths.S2B_MSI_12)]
        metadata['AEC_bands_wl'] = [442.2, 492.3, 559.0, 664.9, 703.8, 739.1, 779.7, 832.9, 864.0, 1610.4, 2185.5]
    elif sensor == 'S2C':
        # source: https://sentinels.copernicus.eu/-/copernicus-sentinel-2c-spectral-response-functions
        metadata['AEC_bands_6S'] = _S2C_RSR()
        metadata['AEC_bands_wl'] = [444.2, 489.0, 560.6, 666.5, 707.1, 741.1, 784.7, 834.6, 865.6, 1612.0, 2191.3]  

    # find central coordinates of MGRS tile 
    m = mgrs.MGRS()
    
    for i, fname in enumerate(files):
        tmp = fname.split('.')
        path = '{}/{}'.format(file,fname)
        
        # granules
        if (fname == 'GRANULE'):
            granules = os.listdir(path)
            
            # check if there is only one granule file 
            n_granule = 0
            
            for granule in granules:
                if granule[0]=='.':continue
                
                n_granule += 1
                if n_granule>1: sys.exit('Warning: more than 1 granule')
                
                metadata['granule'] = '{}/{}/{}/IMG_DATA/'.format(file,fname,granule)
                metadata['MGRS_tile'] = granule.split('_')[1][1:]
                metadata['QI_DATA'] = '{}/{}/{}/QI_DATA'.format(file,fname,granule)
                
                # MGRS
                tile = metadata['MGRS_tile'] + '54905490'
                d = m.toLatLon(tile)
                metadata['lat'] = d[0]
                metadata['lon'] = d[1]
                
                # band files 
                image_files = os.listdir(metadata['granule'])
                for image in image_files: 
                    if image[0]=='.':continue
                    if image[-4:]=='.xml':continue
                    tmp = image.split('_')
                    metadata[tmp[-1][0:3]] = '{}/{}/{}/IMG_DATA/{}'.format(file,fname,granule,image)
                
                # read scene metadata  
                path = '{}/{}/{}/'.format(file,fname,granule)
                granule_files = os.listdir(path)
                for j, grfname in enumerate(granule_files):
                    tmp = grfname.split('.')
                    path = '{}/{}/{}/{}'.format(file,fname,granule,grfname)
                    if (tmp[-1] == ('xml')) & ('MTD' in tmp[0]):
                        xml = tmart.AEC.read_xml_S2(path)
                        metadata.update(xml)
    
    # scaling factors 
    xml2 = tmart.AEC.read_xml_S2_scene(xml_file)
    metadata.update(xml2)
    
    # others 
    metadata['resolution'] = 10
    metadata['reshape_factor'] = int(config['reshape_factor_S2'])
    if metadata['reshape_factor'] % 6 > 0: sys.exit('Warning: reshape_factor_S2 in config must be divisible by 6.')
    
    metadata['window_size'] = int(config['window_size'])
    metadata['AEC_resolution'] = int(config['reshape_factor_S2']) * 10 # resolution 

    metadata['height'] = 10_980
    metadata['width'] = 10_980
    
    # the smallest number divisible by both 6 and reshape_factor
    # 6 comes from 60m/10m resolution, remainder of 6*reshape_f divided by greatest common divisor between 6 and reshape_f
    temp_reshape_factor = (metadata['reshape_factor'] * 6) // math.gcd(metadata['reshape_factor'], 6)
    
    # height and width for AEC, with a few extra rows and columns (ceil = round up)
    metadata['AEC_height'] = math.ceil(metadata['height'] / temp_reshape_factor) * temp_reshape_factor
    metadata['AEC_width'] = math.ceil(metadata['width'] / temp_reshape_factor) * temp_reshape_factor
    
    return metadata 

def _S2C_RSR():
    from Py6S import Wavelength 
    
    # B1
    bands = [Wavelength(0.41, 0.475, [0.0, 0.000881095, 0.00176219, 0.00174064, 0.00171909, 0.0016440399999999998, 0.00156899, 0.001643295, 0.0017176, 0.397548685, 0.79337977, 0.858303195, 0.92322662, 0.9216123, 0.91999798, 0.959226425, 0.99845487, 0.5678938099999999, 0.13733275, 0.069104365, 0.00087598, 0.0007693050000000001, 0.00066263, 0.00064591, 0.00062919, 0.000314595, 0.0]),
             
             # B2
             Wavelength(0.435, 0.545, [0.0, 0.002596715, 0.00519343, 0.004876429999999999, 0.00455943, 0.00436228, 0.00416513, 0.06490035999999999, 0.12563559, 0.514059435, 0.90248328, 0.90103467, 0.89958606, 0.8978790249999999, 0.89617199, 0.8789504, 0.86172881, 0.90781944, 0.95391007, 0.86942361, 0.78493715, 0.812475945, 0.84001474, 0.867130865, 0.89424699, 0.90156546, 0.90888393, 0.941572445, 0.97426096, 0.9114354, 0.84860984, 0.88707073, 0.92553162, 0.911699315, 0.89786701, 0.45806189, 0.01825677, 0.010206815, 0.00215686, 0.00183376, 0.00151066, 0.00147305, 0.00143544, 0.00071772, 0.0]),

             # B3
             Wavelength(0.52, 0.6, [0.0, 0.000416885, 0.00083377, 0.0009703099999999999, 0.00110685, 0.001205135, 0.00130342, 0.00700899, 0.01271456, 0.486878025, 0.96104149, 0.948734295, 0.9364271, 0.942185085, 0.94794307, 0.9735691, 0.99919513, 0.9711927600000001, 0.94319039, 0.9393802600000001, 0.93557013, 0.957859225, 0.98014832, 0.530836185, 0.08152405, 0.04172061, 0.00191717, 0.001490515, 0.00106386, 0.0009255, 0.00078714, 0.00039357, 0.0]),

             # B4
             Wavelength(0.63, 0.705, [0.00072495, 0.000543965, 0.00036298, 0.00042944999999999995, 0.00049592, 0.00149104, 0.00248616, 0.08994199500000001, 0.17739783, 0.5711033400000001, 0.96480885, 0.95213155, 0.93945425, 0.94913285, 0.95881145, 0.97365813, 0.98850481, 0.990942815, 0.99338082, 0.9244100850000001, 0.85543935, 0.444133855, 0.03282836, 0.01725107, 0.00167378, 0.0010631999999999998, 0.00045262, 0.00036716, 0.0002817, 0.00014085, 0.0]),

             # B5
             Wavelength(0.675, 0.74, [0.0, 7.18e-05, 0.0001436, 0.00015978999999999999, 0.00017598, 0.00022334, 0.0002707, 0.00139386, 0.00251702, 0.35224094499999997, 0.70196487, 0.841513695, 0.98106252, 0.98868876, 0.996315, 0.694927965, 0.39354093, 0.1976638, 0.00178667, 0.00102381, 0.00026095, 0.000215055, 0.00016916, 0.00016074, 0.00015232, 7.616e-05, 0.0]),

             # B6
             Wavelength(0.71, 0.77, [0.0, 0.00010817, 0.00021634, 0.000230365, 0.00024439, 0.00043372500000000004, 0.00062306, 0.00705142, 0.01347978, 0.46442181, 0.91536384, 0.95302407, 0.9906843, 0.993840405, 0.99699651, 0.578974745, 0.16095298, 0.081588605, 0.00222423, 0.0013183749999999999, 0.00041252, 0.00033668, 0.00026084, 0.000298975, 0.00033711]),

             # B7
             Wavelength(0.75, 0.82, [0.0, 5.993e-05, 0.00011986, 0.000143435, 0.00016701, 0.00028396, 0.00040091, 0.003271195, 0.00614148, 0.30113213999999994, 0.5961228, 0.7927711, 0.9894194, 0.97679719, 0.96417498, 0.952734625, 0.94129427, 0.710781495, 0.48026872, 0.24303508999999998, 0.00580146, 0.003138325, 0.00047519, 0.000331665, 0.00018814, 0.00016361000000000001, 0.00013908, 6.954e-05, 0.0]),

             # B8
             Wavelength(0.76, 0.925, [0.0, 0.000624345, 0.00124869, 0.0019381200000000002, 0.00262755, 0.00578735, 0.00894715, 0.030729779999999998, 0.05251241, 0.225869775, 0.39922714, 0.6761808949999999, 0.95313465, 0.9755721900000001, 0.99800973, 0.984465605, 0.97092148, 0.94309009, 0.9152587, 0.8840158149999999, 0.85277293, 0.8236130500000001, 0.79445317, 0.770022035, 0.7455909, 0.7287530499999999, 0.7119152, 0.7047625449999999, 0.69760989, 0.6953143799999999, 0.69301887, 0.69806925, 0.70311963, 0.698928445, 0.69473726, 0.6777183950000001, 0.66069953, 0.6447894949999999, 0.62887946, 0.61626223, 0.603645, 0.5860471549999999, 0.56844931, 0.5543814149999999, 0.54031352, 0.532305685, 0.52429785, 0.515122295, 0.50594674, 0.49497073999999996, 0.48399474, 0.465841215, 0.44768769, 0.43347445, 0.41926121, 0.33062718, 0.24199315, 0.14207378999999998, 0.04215443, 0.025086209999999998, 0.00801799, 0.0051807549999999996, 0.00234352, 0.0016449799999999999, 0.00094644, 0.00047322, 0.0]),

             # B8A
             Wavelength(0.83, 0.9, [0.0, 8.291e-05, 0.00016582, 0.00018964, 0.00021346, 0.00032219, 0.00043092, 0.002265495, 0.00410007, 0.15940945, 0.31471883, 0.6560948, 0.99747077, 0.980849005, 0.96422724, 0.950272335, 0.93631743, 0.7946804599999999, 0.65304349, 0.33604138499999997, 0.01903928, 0.010021164999999999, 0.00100305, 0.0006311349999999999, 0.00025922, 0.000206555, 0.00015389, 7.6945e-05, 0.0]),

             # B11
             Wavelength(1.54, 1.68, [0.00063023, 0.0018798299999999999, 0.00312943, 0.005263185, 0.00739694, 0.01565511, 0.02391328, 0.06117299, 0.0984327, 0.23930785999999998, 0.38018302, 0.59450675, 0.80883048, 0.874626765, 0.94042305, 0.943369175, 0.9463153, 0.94972684, 0.95313838, 0.9564766899999999, 0.959815, 0.962702305, 0.96558961, 0.96870959, 0.97182957, 0.972485915, 0.97314226, 0.9739119650000001, 0.97468167, 0.976320635, 0.9779596, 0.97931323, 0.98066686, 0.97925894, 0.97785102, 0.97766712, 0.97748322, 0.98082709, 0.98417096, 0.985082885, 0.98599481, 0.98719095, 0.98838709, 0.988893255, 0.98939942, 0.832594345, 0.67578927, 0.44498283, 0.21417639, 0.140279975, 0.06638356, 0.053024379999999996, 0.0396652, 0.03049049, 0.02131578, 0.01065789, 0.0]),

             # B12
             Wavelength(2.065, 2.315, [0.0, 0.023090065, 0.04618013, 0.048921725, 0.05166332, 0.05459978, 0.05753624, 0.065419965, 0.07330369, 0.092886035, 0.11246838, 0.167354715, 0.22224105, 0.34471193499999997, 0.46718282, 0.62386936, 0.7805559, 0.8556526849999999, 0.93074947, 0.944066885, 0.9573843, 0.95819447, 0.95900464, 0.95755843, 0.95611222, 0.95489043, 0.95366864, 0.95641225, 0.95915586, 0.96746219, 0.97576852, 0.9800847850000001, 0.98440105, 0.977876725, 0.9713524, 0.96981841, 0.96828442, 0.980440155, 0.99259589, 0.982571395, 0.9725469, 0.93742795, 0.902309, 0.8955609250000001, 0.88881285, 0.88684352, 0.88487419, 0.88250849, 0.88014279, 0.873219025, 0.86629526, 0.85450101, 0.84270676, 0.824489185, 0.80627161, 0.80031187, 0.79435213, 0.80325144, 0.81215075, 0.83221284, 0.85227493, 0.8825403949999999, 0.91280586, 0.92192082, 0.93103578, 0.93758567, 0.94413556, 0.94434701, 0.94455846, 0.943006045, 0.94145363, 0.94389716, 0.94634069, 0.94792796, 0.94951523, 0.9488188200000001, 0.94812241, 0.9491464249999999, 0.95017044, 0.9554948999999999, 0.96081936, 0.9699252899999999, 0.97903122, 0.959971275, 0.94091133, 0.8396976300000001, 0.73848393, 0.577217825, 0.41595172, 0.30558528, 0.19521884, 0.14127084, 0.08732284, 0.065392885, 0.04346293, 0.033471429999999996, 0.02347993, 0.0183637, 0.01324747, 0.006623735, 0.0])]
    
    return bands 

