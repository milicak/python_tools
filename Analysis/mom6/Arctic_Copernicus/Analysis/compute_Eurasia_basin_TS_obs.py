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
df = xr.open_dataset(root_folder + 'TS_woa_Arctic_Copernicus.nc',decode_times=False)
gridname = root_folder + 'ocean_geometry.nc'
gr = xr.open_dataset(gridname)

# regions = loadmat('../../../mitgcm/Arctic_4km/Analysis/region_masks.mat')
# vertices = np.transpose(np.array([regions['lon7'].flatten(),regions['lat7'].flatten()]))
lon1 = np.loadtxt('Eurasia_lon.txt')
lat1 = np.loadtxt('Eurasia_lat.txt')
vertices = np.transpose(np.array([lon1,lat1]))
mpath = Path( vertices ) # the vertices of the polygon

lontemp = (gr.geolon+180)%360-180
lonlat = np.dstack((lontemp,gr.geolat))
lonlat_flat = lonlat.reshape((-1, 2))
mask_flat = mpath.contains_points(lonlat_flat)
mask = mask_flat.reshape(gr.geolon.shape)
mask = mask*1
mask[785:830,1175:1230] = 1

depthmask = xr.where(gr.D*mask>500,1,0)
depthmask = depthmask.rename({'lonh':'xh','lath':'yh'})
area = gr.Ah.rename({'lonh':'xh','lath':'yh'})
df = df.rename({'lonh':'xh','lath':'yh'})

tmp = area*depthmask
tmp = tmp.fillna(0)
dnm = np.tile(tmp,(102,1,1))
dnm[np.isnan(df.temp_woa)]=np.nan
ds = df.temp_woa*dnm
ds = ds.fillna(0)
temp1 = np.zeros(102)
for kind in range(0,102):
    temp1[kind] = ds[kind,:,:].sum()/np.nansum(dnm[kind,:,:])


ds = df.salt_woa*dnm
ds = ds.fillna(0)
salt1 = np.zeros(102)
for kind in range(0,102):
    salt1[kind] = ds[kind,:,:].sum()/np.nansum(dnm[kind,:,:])

# create dataset
dfs = xr.Dataset({
    'salt_Eurasia': xr.DataArray(
                data   = salt1,
                dims   = ['depth'],
                coords = {'depth': np.copy(ds.depth)},
                attrs  = {
                    'units'     : 'psu'
                    }
                ),
    'temp_Eurasia': xr.DataArray(
                data   = temp1,
                dims   = ['depth'],
                coords = {'depth': np.copy(ds.depth)},
                attrs  = {
                    'units'     : 'C'
                    }
                )
            },
    )

dfs.to_netcdf('Eurasia_TS_WOA13.nc')
