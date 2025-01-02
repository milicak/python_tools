import os
import numpy as np
import xesmf as xe
import xarray as xr
# import scipy.io
# from scipy.io import savemat
# from scipy.io import loadmat
# import matplotlib.colors as colors
# from scipy.signal import medfilt2d
# import netCDF4
# import matplotlib.pyplot as plt
# from scipy.interpolate import griddata
# from matplotlib.path import Path
#for interpolation
# from scipy.spatial import cKDTree
from HCtFlood.kara import flood_kara
# from PyCNAL_regridding import *


root_folder = '/okyanus/users/milicak//dataset/NEMO/ArabianSea/'
ls1 = sorted(glob.glob('/okyanus/users/milicak/dataset/obs/woa13_decav_s*.nc'))
ls1 = ls1[1:-1]
dds = xr.open_mfdataset(ls1,decode_times=False)
path_regional_grid = root_folder + './domain_cfg.nc'

variable = 's_an'
ds = dds[variable][:,0,:,:]
ds = ds.to_dataset(name=variable)
df = flood_kara(ds[variable][:,:,:], xdim='lon', ydim='lat', zdim='depth')
df = df.to_dataset(name='sss')

df2 = xr.open_dataset(path_regional_grid)
variables=['nav_lon','nav_lat']
ds2 = df2[variables]
ds2 = ds2.rename({'nav_lon':'lon','nav_lat':'lat'})

df.load()

# build regridder
regridder = xe.Regridder(df, ds2, 'nearest_s2d')

#apply regridder
dr_out = regridder(df['sss'])
# dr_out = dr_out.ffill('depth')
tmp = dr_out.to_dataset(name='sss')


# vertical 1D interpolation
# tmp = tmp.drop('time_counter')
tmp = tmp.rename({'time': 'time_counter'})
all_vars = list(tmp.data_vars.keys()) + list(tmp.coords.keys())
encodings = {v: {'_FillValue': None} for v in all_vars}
encodings['time_counter'].update({'dtype':'float64', 'calendar': 'gregorian'})
ftmp = root_folder + 'SSS_climatology_WOA18.nc'
tmp.to_netcdf(
        ftmp,
        format='NETCDF4_CLASSIC',
        engine='netcdf4',
        unlimited_dims=['time_counter'])
tmp.close()
