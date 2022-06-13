import os
import numpy as np
import xesmf as xe
import xarray as xr
import scipy.io
from scipy.io import savemat
from scipy.io import loadmat
import matplotlib.colors as colors
from scipy.signal import medfilt2d
import netCDF4
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.path import Path
#for interpolation
from scipy.spatial import cKDTree
from HCtFlood.kara import flood_kara
from os.path import exists

mom_dir = '/okyanus/users/milicak/dataset/MOM6/NA12/'
path_regional_grid = mom_dir + './ocean_hgrid.nc'
df2 = xr.open_dataset(path_regional_grid)
lon_rho = np.copy(df2['x'][1::2,1::2])
lat_rho = np.copy(df2['y'][1::2,1::2])
nj,ni = lon_rho.shape
ds2 = df2['x'][1::2,1::2]
ds2 = ds2.to_dataset(name='lon')
ds2['lat']=df2['y'][1::2,1::2]
ds2 = ds2.rename_dims({'nxp': 'x','nyp': 'y'})

root_folder = '/okyanus/users/milicak/dataset/SODA3_12_21/monthlyfiles/'
f = 'soda3.12.2_mn_ocean_reg'

year = 1996
current = xr.open_dataset(f"{root_folder}/{f}_{year}.nc")
next_data = xr.open_dataset(f"{root_folder}/{f}_{year+1}.nc").isel(time=0)
previous = xr.open_dataset(f"{root_folder}/{f}_{year-1}.nc").isel(time=-1)
out = xr.concat([previous, current, next_data], dim="time")
out = out.rename({'xt_ocean': 'lon', 'yt_ocean': 'lat','st_ocean':'depth'})

# swap dimensions
out2 = out.temp.transpose("time","depth","lat","lon")
out2 = out2.to_dataset(name='temp')
out2['salinity'] = out.salt.transpose("time","depth","lat","lon")
out2['ssh'] = out.ssh.transpose("time","lat","lon")

# build regridder
regridder2 = xe.Regridder(out2, ds2, 'patch',reuse_weights=False, periodic=True)


