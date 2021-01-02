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

regions = loadmat('region_masks.mat')
vertices = np.transpose(np.array([regions['lon1'].flatten(),regions['lat1'].flatten()]))
mpath = Path( vertices ) # the vertices of the polygon

root_folder = '/archive/milicak/MITgcm_c65/Projects/Arctic_4km/'

expid = 'Exp02_1';

gridname = root_folder + expid + '/grid.nc'
gr = xr.open_dataset(gridname)

lonlat = np.dstack((gr.XC,gr.YC))
lonlat_flat = lonlat.reshape((-1, 2))
mask_flat = mpath.contains_points(lonlat_flat)
mask = mask_flat.reshape(gr.XC.shape)
mask = mask*1

fyear = 1992
# fyear = 1992
lyear = 2018
# lyear = fyear+1
datadir = root_folder+expid
# os.chdir(datadir)
prename = '2DArcticOcean_'
postname = '_avg.nc'

fname = datadir+'/'+prename+'*'+'SIarea*'
list=sorted(glob.glob(fname))

time = pd.date_range('1992-01-01', freq='M', periods=12 * 25)

df = xr.open_mfdataset(list)
# df['time'] = time

si = df*gr.rA*mask
SI = si.sum(dim=['i','j'])
fname = root_folder + 'ncfiles/' + expid + '_seaice_area_Barents.nc'
SI.to_netcdf(fname)
siann = SI.groupby('time.year').mean('time')*1e-12

aa = xr.where(df.SIarea<0.15,0,1)
aa = aa*gr.rA*mask
SIext = aa.sum(dim=['i','j'])
ds = SIext.to_dataset(name='SI_extent')
fname = root_folder + 'ncfiles/' + expid + '_seaice_extent_Barents.nc'
ds.to_netcdf(fname)

