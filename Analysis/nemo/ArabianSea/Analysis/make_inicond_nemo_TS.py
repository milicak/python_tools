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
fname1 = 'cmems_mod_glo_phy_my_0.083_P1D-m_1699373113831.nc'
dds = xr.open_dataset(root_folder + fname1)
path_regional_grid = root_folder + './domain_cfg.nc'

variable = 'thetao'
ds = dds[variable][0,:,:,:]
ds = ds.to_dataset(name='thetao')
df = flood_kara(ds[variable][:,:,:], xdim='longitude', ydim='latitude', zdim='depth')
df = df.to_dataset(name='temp')
df = df.rename({'longitude':'lon','latitude':'lat'})

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
ds = dds[variable][0,:,:,:]
ds = ds.to_dataset(name=variable)
df = flood_kara(ds[variable][:,:,:], xdim='longitude', ydim='latitude', zdim='depth')
df = df.to_dataset(name=variable)
df = df.rename({'longitude':'lon','latitude':'lat'})
df.load()
#apply regridder
dr_out = regridder(df[variable])
dr_out = dr_out.ffill('depth')
df_out[variable] = dr_out

# vertical 1D interpolation
zr = df2.e3t_1d[0,:].cumsum()
tmp = df_out.interp(depth=zr)
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
