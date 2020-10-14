import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
# import ESMF
# from mpl_toolkits.basemap import Basemap                                            
import geopy.distance
#import cartopy.crs as ccrs                                                          
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
 

ds = xr.open_dataset('uTSS_lobc_chunk_0640.nos.nc')
df = pd.read_excel('20190220_1033_ODV_TMAMVerisi.xlsx',sheet_name='2017data')
lons = np.copy(ds.longitude)
lats = np.copy(ds.latitude)
utss_tree = cKDTree(zip(lons,lats))

timecol = 'yyyy-mm-ddThh:mm:ss.sss'
depthcol = 'Depth [m]'
tempcol = 'Temperature [degC]'
saltcol = 'Salinity [psu]'

stations = df[timecol].unique()

grade = 4

# for each station
timeind = 0
# expname = '2017-01_ISKI_MARMARA'
expname = '2017-08_BUTUNLESIK_MARMARA'
tmp_pro = df.loc[df['Cruise']==expname]
# stationame = 'M23'  # 'M23' 'B2' 'B2B' 'B5' 'B7' 'K0A' 'K0' 'M14' 'M18' 'K0' 'K2' 'BL1' 'B14' 'B13' 'M1' 'M8' 'MBA' 'M20'
# stationames = ['M23', 'B2', 'B2B', 'B5', 'B7', 'K0A', 'K0', 'M14']
# stationames = ['M20', 'M23', 'M8', 'M14', 'MY2', 'MBA', 'M1', 'SB', 'B7',
#                'B13', 'B14']
stationames = ['K0', '45C', 'MD102', 'M20', 'MD104', 'M14A', 'YSA', 'DIPTAR1Y',
               'DIPTAR2Y', 'DIPTAR3Y', 'DIPTAR4Y', 'MD75', 'AR1', 'MD18',
               'MD101',  'MD13A', 'KD1', 'MD67', 'MD10A', 'MD10B']

for stationame in stationames:
    print(stationame)
    tmp_prof = tmp_pro.loc[tmp_pro['Station']==stationame]
    lon_thlweg = np.array([ tmp_prof['Longitude'].tolist()[0]])
    lat_thlweg = np.array([ tmp_prof['Latitude'].tolist()[0]])
    # coordinates of K0 (lat, lon)
    coords_tmp = np.concatenate((lat_thlweg, lon_thlweg))
    d_utss, inds_utss = utss_tree.query([lon_thlweg[0],lat_thlweg[0]],k = grade)
    df2 = pd.DataFrame({"d_utss": d_utss, "inds_utss": inds_utss})
    fname = 'stations_indices/'+stationame+'_obs.csv'
    df2.to_csv(fname)


