# This file is part of TMart.
#
# Copyright 2024 Yulun Wu.
#
# TMart is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.


# Download hourly ancillary data from NASA Ocean Color 
def anci_download(file, basefile, username, password):
    
    import requests, time, os, sys
    url = 'https://oceandata.sci.gsfc.nasa.gov/ob/getfile/' + basefile
    auth = (username, password)
    local_file = os.path.join(file, 'tmart_ancillary', basefile)
    
    if os.path.exists(local_file):
        print('File exists: ' + local_file)
        pass
    
    else:
        print('Downloading: ' + local_file)
        try:     
            if os.path.exists(os.path.dirname(local_file)) is False:
                os.makedirs(os.path.dirname(local_file))
            
            # source: https://github.com/acolite/acolite/blob/main/acolite/shared/download_file.py
            with requests.Session() as session:
                r1 = session.request('get', url, verify=True)
                r = session.get(r1.url, auth=auth, verify=True)
                time.sleep(1)
            
                if (r.ok):
                    with open(local_file, 'wb') as f:
                        for chunk in r.iter_content(chunk_size=1024*1024):
                            if chunk: f.write(chunk)
                  
        except: 
            sys.exit('Warning: failed to download: ' + str(basefile) + '. Please check your credentials and if file is available on NASA Ocean Color website. Note that the files are typically available a month after image acquisition.')

    if os.path.exists(local_file):
        return local_file       
    else:
        sys.exit('Warning: failed to download: ' + str(basefile) + '. Please check your credentials and if file is available on NASA Ocean Color website. Note that the files are typically available a month after image acquisition.')

