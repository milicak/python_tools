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


root_folder2 = '/okyanus/users/milicak//dataset/NEMO/TSS_BS_AGRIF/'
fname1 = 'cmems_1993_01_01_T.nc'
dds = xr.open_dataset(root_folder2 + fname1)
root_folder = '/okyanus/users/milicak//dataset/NEMO/TSS_AGRIF/'
path_regional_grid = root_folder + './domain_cfg.nc'

variable = 'thetao'
ds = dds[variable][0,:,:,:]
ds = ds.to_dataset(name='thetao')
ds = ds.ffill('depth')
df = flood_kara(ds[variable][:80,:,:], xdim='lon', ydim='lat', zdim='depth')
df = df.to_dataset(name='temp')

df2 = xr.open_dataset(path_regional_grid)
variables=['nav_lon','nav_lat']
ds2 = df2[variables]
ds2 = ds2.rename({'nav_lon':'lon','nav_lat':'lat'})

df.load()

# build regridder
regridder = xe.Regridder(df, ds2, 'nearest_s2d')

#apply regridder
dr_out = regridder(df['temp'])
dr_out = dr_out.ffill('depth')
df_out = dr_out.to_dataset(name='thetao')

# salt
variable = 'so'
fname1 = 'cmems_1993_01_01_S.nc'
dds = xr.open_dataset(root_folder2 + fname1)
ds = dds[variable][0,:,:,:]
ds = ds.to_dataset(name=variable)
df = flood_kara(ds[variable][:80,:,:], xdim='lon', ydim='lat', zdim='depth')
df = df.to_dataset(name=variable)
df.load()
#apply regridder
dr_out = regridder(df[variable])
dr_out = dr_out.ffill('depth')
df_out[variable] = dr_out

# vertical 1D interpolation
zr = df2.e3t_1d[0,:].cumsum()
tmp = df_out.interp(depth=zr)
tmp = tmp.ffill('nav_lev')
tmp = tmp.bfill('nav_lev')
tmp = tmp.drop('time_counter')
tmp = tmp.rename({'time': 'time_counter'})
all_vars = list(tmp.data_vars.keys()) + list(tmp.coords.keys())
encodings = {v: {'_FillValue': None} for v in all_vars}
encodings['time_counter'].update({'dtype':'float64', 'calendar': 'gregorian'})
ftmp = root_folder + 'TSy1993_m01_d01.nc'
tmp.to_netcdf(
        ftmp,
        format='NETCDF4_CLASSIC',
        engine='netcdf4',
        unlimited_dims=['time_counter'])
tmp.close()
