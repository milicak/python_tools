import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
# import ESMF
# from mpl_toolkits.basemap import Basemap                                            
import geopy.distance
#import cartopy.crs as ccrs                                                          
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
 
plt.ion()




df = xr.open_dataset('/work/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/OUT/uTSS_lobc_chunk_0700.ous.nc')
# coordinates of K0 (lat, lon)
coords_K0 = (41.226, 29.135)
# coordinates of K2 (lat, lon)
coords_K2 = (41.285, 29.180)
# coordinates of M23 (lat, lon)
coords_M23 = (40.70, 28.78348)

dst = 1e6*np.ones(df.latitude.shape)
for ind in range(0,np.copy(df.longitude.shape)-1):
    coords_1=(df.latitude[ind],df.longitude[ind])
    dst[ind] = geopy.distance.geodesic(coords_1, coords_M23).km
    # dst[ind] = geopy.distance.geodesic(coords_1, coords_K0).km
    # dst[ind] = geopy.distance.geodesic(coords_1, coords_K2).km
    # dst[ind] = geopy.distance.vincenty(coords_1, coords_K0).km


index_min = np.argmin(dst)
print(index_min)


# dst = 1e6*np.ones(df.latitude.shape)
# for ind in range(0,np.copy(df.longitude.shape)-1):
#     coords_1=(df.latitude[ind],df.longitude[ind])
#     dst[ind] = geopy.distance.geodesic(coords_1, coords_K2).km
#
#
# index_min = np.argmin(dst)
# print(index_min)
