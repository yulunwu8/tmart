# This file is part of TMart.
#
# Copyright 2024 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


def anci_list_files(metadata):
    from datetime import datetime, timedelta
    
    date_string = metadata['time'][0:13]
    date_format = '%Y-%m-%dT%H'
    
    sat_overpass = datetime.strptime(date_string, date_format)
    sat_overpass_h0 = sat_overpass.replace(second=0, microsecond=0, minute=0)
    sat_overpass_h1 = sat_overpass_h0 + timedelta(hours=1)
    
    # isodate = sat_overpass.strftime("%Y-%m-%d")
    # year = sat_overpass.strftime("%Y")
    # jday = sat_overpass.strftime("%j").zfill(3)
    # yjd = "{}{}".format(sat_overpass.strftime("%Y"),jday)
    
    # dtime_next = sat_overpass + timedelta(days=1)
    # year_next = dtime_next.strftime("%Y")
    # jday_next = dtime_next.strftime("%j").zfill(3)
    # yjd_next = "{}{}".format(year_next,jday_next)
    
    # file_ozone = "N{}00_O3_AURAOMI_24h.hdf".format(yjd)
    # test = ["N{}{}_MET_NCEP_6h.hdf".format(yjd,h) for h in ['00','06','12','18']]
    # test = "N{}00_MET_NCEP_6h.hdf".format(yjd_next)
    
    # Aerosol 
    aer_h0 = 'GMAO_MERRA2.' + sat_overpass_h0.strftime("%Y%m%dT%H") + '0000.AER.nc'
    aer_h1 = 'GMAO_MERRA2.' + sat_overpass_h1.strftime("%Y%m%dT%H") + '0000.AER.nc'
    
    # Ozone and water vapour 
    met_h0 = 'GMAO_MERRA2.' + sat_overpass_h0.strftime("%Y%m%dT%H") + '0000.MET.nc'
    met_h1 = 'GMAO_MERRA2.' + sat_overpass_h1.strftime("%Y%m%dT%H") + '0000.MET.nc'   
    
    return {'AER': [aer_h0,aer_h1],'MET': [met_h0,met_h1]}
    
