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


root_folder = '/okyanus/users/milicak//dataset/NEMO/TSS_BS_AGRIF/'
fname1 = 'T_utss-sdc2019m01-sol.nc'
dds = xr.open_dataset(root_folder + fname1)
path_regional_grid = root_folder + './domain_cfg.nc'

gr = xr.open_dataset('/okyanus/users/milicak//dataset/NEMO/BS_AGRIF/domain_cfg_noagrif.nc')
mask = xr.where(gr.bathy_metry!=0,1,0);
mask = mask.isel(x=slice(0,100),y=slice(0,26))
mask = mask.rename({'x':'X','y':'Y'})

variable = 'votemper'
# df = dds[variable][0,:,:,:]
# ds = dds[variable][0,:,:,:]
ds = dds.isel(X=slice(0,100),Y=slice(0,26))
mask['time_counter']  = ds['time_counter']
ds.load()
ds = ds.where(mask!=0)
# ds = ds.to_dataset(name='thetao')
df = flood_kara(ds[variable][:,:,:], xdim='X', ydim='Y', zdim='Z')
df = df.to_dataset(name='temp')
# df = dds.isel(X=slice(0,100),Y=slice(0,26))
# df = df.rename({'X':'lon','Y':'lat'})

df2 = xr.open_dataset(path_regional_grid)
variables=['nav_lon','nav_lat']
ds2 = df2[variables]
ds2 = ds2.rename({'nav_lon':'lon','nav_lat':'lat'})

df.load()
df = df.rename({'X':'lon','Y':'lat'})

# build regridder
regridder = xe.Regridder(df, ds2, 'nearest_s2d')

#apply regridder
dr_out1 = regridder(df['temp'])
dr_out1 = dr_out1.ffill('Z')
df_out1 = dr_out1.to_dataset(name='thetao')

# salt
variable = 'vosaline'
fname1 = 'S_utss-sdc2019m01-sol.nc'
dds = xr.open_dataset(root_folder + fname1)
ds = dds.isel(X=slice(0,100),Y=slice(0,26))
ds.load()
ds = ds.where(mask!=0)
# ds = ds.to_dataset(name='thetao')
df = flood_kara(ds[variable][:,:,:], xdim='X', ydim='Y', zdim='Z')
df = df.to_dataset(name='salt')
df.load()
df = df.rename({'X':'lon','Y':'lat'})

#apply regridder
dr_out = regridder(df['salt'])
dr_out = dr_out.ffill('Z')
df_out1[variable] = dr_out
# over write for some values in deep
for k in range(110,121):
    df_out1.thetao[0,k,:,:] = df_out1.thetao[0,109,:,:]
    df_out1.vosaline[0,k,:,:] = df_out1.vosaline[0,109,:,:]

# vertical 1D interpolation
zr = df2.e3t_1d[0,:].cumsum()
zr = zr.drop('time_counter')
gr1 = xr.open_dataset('~/dataset/NEMO/BS/domain_cfg.nc')
zr_bs = gr1.e3t_1d[0,:].cumsum()
df_out1['Z'] = zr_bs
tmp = df_out1.interp(Z=zr)
tmp = tmp.ffill('nav_lev')
tmp = tmp.rename({'time': 'time_counter'})
tmp['time_counter']  = ds['time_counter']

all_vars = list(tmp.data_vars.keys()) + list(tmp.coords.keys())
encodings = {v: {'_FillValue': None} for v in all_vars}
encodings['time_counter'].update({'dtype':'float64', 'calendar': 'gregorian'})
ftmp = root_folder + 'TS_climatology_Marmara.nc'
tmp.to_netcdf(
        ftmp,
        format='NETCDF4_CLASSIC',
        engine='netcdf4',
        unlimited_dims=['time_counter'])
tmp.close()

