import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import os
from datetime import date
from matplotlib.path import Path
from scipy.io import loadmat
from pandas.tseries.offsets import DateOffset

# lon1,lat1 is for Kara and Barents Sea
# lon2,lat2 is for Greenland Sea
# lon3,lat3 is for Hudson Bay
# lon4,lat4 is for CAA
# lon5,lat5 is for Arctic Ocean Canadian side
# lon6,lat6 is for Labrador Sea/ Baffin Bay
# lon7,lat7 is for Arctic Ocean Eurasian side
# lon8,lat8 is for Bering Sea
# lon9,lat9 is for Chukchi Sea
# lon10,lat10 is for East Siberian Sea

root_folder = '/archive/milicak/MOM6-examples/Projects/Arctic_Copernicus/'
gridname = root_folder + 'ocean_geometry.nc'
gr = xr.open_dataset(gridname)

# regions = loadmat('../../../mitgcm/Arctic_4km/Analysis/region_masks.mat')
# vertices = np.transpose(np.array([regions['lon7'].flatten(),regions['lat7'].flatten()]))
lon1 = np.loadtxt('Canada_lon.txt')
lat1 = np.loadtxt('Canada_lat.txt')
aa = lon1[lon1<0]+360
xx = np.append(lon1[0:4],aa)
xx = np.append(xx,lon1[-6:])
vertices = np.transpose(np.array([xx,lat1]))
mpath = Path( vertices ) # the vertices of the polygon

lontemp = gr.geolon
lonlat = np.dstack((lontemp,gr.geolat))
lonlat_flat = lonlat.reshape((-1, 2))
mask_flat = mpath.contains_points(lonlat_flat)
mask = mask_flat.reshape(gr.geolon.shape)
mask = mask*1

depthmask = xr.where(gr.D*mask>500,1,0)
depthmask = depthmask.rename({'lonh':'xh','lath':'yh'})

ls1 = sorted(glob.glob(root_folder + '*ocean_month*'))
ls2 = ls1[12:24]
for itr in range(3,240,2):
    ls2.extend(ls1[itr*12:itr*12+12])

df = xr.open_mfdataset(ls2[0],decode_times=False)
area = gr.Ah.rename({'lonh':'xh','lath':'yh'})
tmp = area*depthmask
tmp = tmp.fillna(0)
dnm = np.tile(tmp,(df.thetao.shape[1],1,1))
dnm[np.isnan(df.thetao[0,:,:,:])]=np.nan
temp1 = np.zeros((len(ls2),df.thetao.shape[1]))
for tind in range(0,len(ls2)):
    print(tind)
    df = xr.open_dataset(ls2[tind],decode_times=False)
    ds = df.thetao*dnm
    ds = ds.fillna(0)
    for kind in range(0,df.thetao.shape[1]):
        temp1[tind,kind] = ds[0,kind,:,:].sum()/np.nansum(dnm[kind,:,:])


salt1 = np.zeros((len(ls2),df.so.shape[1]))
for tind in range(0,len(ls2)):
    print(tind)
    df = xr.open_dataset(ls2[tind],decode_times=False)
    ds = df.so*dnm
    ds = ds.fillna(0)
    for kind in range(0,df.so.shape[1]):
        salt1[tind,kind] = ds[0,kind,:,:].sum()/np.nansum(dnm[kind,:,:])


# set time
time = pd.date_range("1996-01-15", freq=DateOffset(months=1), periods=temp1.shape[0])
dfs = xr.Dataset({
    'salt_Canada': xr.DataArray(
                data   = salt1,
                dims   = ['time','depth'],
        coords = {'time': time, 'depth': np.copy(ds.z_l)},
                attrs  = {
                    'units'     : 'psu'
                    }
                ),
    'temp_Canada': xr.DataArray(
                data   = temp1,
                dims   = ['time','depth'],
        coords = {'time': time, 'depth': np.copy(ds.z_l)},
                attrs  = {
                    'units'     : 'C'
                    }
                )
            },
    )

dfs.to_netcdf('Canada_TS_mom6.nc')

