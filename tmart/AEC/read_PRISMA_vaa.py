# This file is part of T-Mart.
#
# Copyright 2024 Yulun Wu.
#
# T-Mart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# Read PRSIMA vaa

def read_PRISMA_vaa(file):
    
    import h5py
    import numpy as np

    # From Dr. Claudia Giardino, giardino.c@irea.cnr.it
    def view_angles(xPos, yPos, zPos, lon , lat , h=0):
        phiRad = lat * np.pi / 180;
        lamRad = lon * np.pi / 180;

        a_ax = 6378137.
        b_ax = 6356752.314205;
        e2 = 1 - (np.power(b_ax, 2) / np.power(a_ax, 2));
        N = a_ax / np.sqrt(1 - (e2 * np.power(np.sin(phiRad), 2)));

        #/*From WGS84 to ECEF*/
        x = (N + h) * np.cos(phiRad) * np.cos(lamRad);
        y = (N + h) * np.cos(phiRad) * np.sin(lamRad);
        z = (N*(1 - e2) + h)*np.sin(phiRad);
        
        #/*Computing normal vector*/
        n = np.array([np.cos(phiRad)*np.cos(lamRad), np.cos(phiRad)*np.sin(lamRad), np.sin(phiRad)]);
        
        #/*Computing east vector*/
        E = np.array([-np.sin(lamRad), np.cos(lamRad), 0]);
        
        #/*Computing the north vector*/
        N1 = n[1]*E[2] - n[2]*E[1];
        N2 = n[2]*E[0] - n[0]*E[2];
        N3 = n[0]*E[1] - n[1]*E[0];
        North = np.array([N1, N2, N3]);
        
        #/*Generate the satellite range*/
        v_sc = np.array([ xPos-x, yPos-y, zPos-z ]);
        r_sc = np.sqrt(np.power(v_sc[0], 2) + np.power(v_sc[1], 2) + np.power(v_sc[2], 2));
        v_norm =np.array([ 1/r_sc* v_sc[0], 1/r_sc * v_sc[1], 1/r_sc * v_sc[2]]);

        #/*ZENITH*/
        tmpZ = n[0]*v_norm[0] + n[1]*v_norm[1] + n[2]*v_norm[2];
        zenith = np.arccos(tmpZ) * 180/np.pi;

        #/*AZIMUTH*/
        I = v_sc[0]*E[0] + v_sc[1]*E[1] + v_sc[2]*E[2];
        m = v_sc[0]*North[0] + v_sc[1]*North[1] + v_sc[2]*North[2];
        azimuth = np.arctan2(I,m)*180/np.pi; 
        
        return zenith, azimuth

    with h5py.File(file, mode='r') as f:
        
        a_pos_x = f['Info']['Ancillary']['PVSdata']['Wgs84_pos_x'][:]
        a_pos_y = f['Info']['Ancillary']['PVSdata']['Wgs84_pos_y'][:]
        a_pos_z = f['Info']['Ancillary']['PVSdata']['Wgs84_pos_z'][:]
    
        src = 'HCO'
    
        lat_key = 'Latitude_SWIR'
        lon_key = 'Longitude_SWIR'
    
        a_lat = f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Geolocation Fields'][lat_key][:]
        a_lon = f['HDFEOS']['SWATHS']['PRS_L1_{}'.format(src)]['Geolocation Fields'][lon_key][:]
        
        pos_x = np.nanmean(a_pos_x)
        pos_y = np.nanmean(a_pos_y)
        pos_z = np.nanmean(a_pos_z)
        
        lat = np.nanmean(a_lat)
        lon = np.nanmean(a_lon)
        
    angle = view_angles(pos_x, pos_y, pos_z, lon , lat)
    
    return angle

