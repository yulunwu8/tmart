
# Test adjacency effect correction 

import tmart

file = '/Volumes/San/test/script/S2B_MSIL1C_20200812T000249_N0209_R030_T55HFV_20200812T011745.SAFE'


# NASA EarthData Credentials, you may need to approve OB.DAAC Data Access in your EarthData account.
username = 'abcdef'
password = '123456'


if __name__ == "__main__":
    test = tmart.AEC.run(file, username, password, overwrite=True, n_photon = 10_000)



