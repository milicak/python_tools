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


root_folder = '/okyanus/users/milicak/dataset/SODA3_12_21/inicon_sodafiles/'
fname1 = 'soda3.12.2_5dy_ocean_reg_1980_01_03.nc'
dds = xr.open_dataset(root_folder + fname1)
mom_dir = '/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/'
path_regional_grid = mom_dir + './ocean_hgrid.nc'

variable = 'temp'

# I set only high lat values so that in the deep ocean flood would work
ds = dds[variable][0,:,210:,:]
ds = ds.to_dataset(name=variable)
df = flood_kara(ds[variable][:,:,:], xdim='xt_ocean', ydim='yt_ocean',
                zdim='st_ocean')
df = df.to_dataset(name='temp')
df = df.rename({'xt_ocean': 'lon', 'yt_ocean': 'lat'})

df2 = xr.open_dataset(path_regional_grid)
lon_rho = np.copy(df2['x'][1::2,1::2])
lat_rho = np.copy(df2['y'][1::2,1::2])
nj,ni = lon_rho.shape
ds2 = df2['x'][1::2,1::2]
ds2 = ds2.to_dataset(name='lon')
ds2['lat']=df2['y'][1::2,1::2]
ds2 = ds2.rename_dims({'nxp': 'x','nyp': 'y'})

# for regridding
df = df.rename({'st_ocean':'depth'})
temp = np.copy(df.temp)
dstemp = xr.DataArray(temp, coords=[df.time, df.depth, df.lat, df.lon],
                            dims=["time","depth","lat","lon"])
dstemp = dstemp.ffill(dim='depth')
df['temp'] = dstemp

# build regridder
regridder = xe.Regridder(df, ds2, 'nearest_s2d')

#apply regridder
dr_out = regridder(df['temp'])
# create 3d ssh variable
# not sure if i need to multiply with mask
var = np.copy(dr_out)
nk = var.shape[1]


# salinity
variable = 'salt'
ds = dds[variable][0,:,210:,:]
ds = ds.to_dataset(name=variable)
df = flood_kara(ds[variable][:,:,:], xdim='xt_ocean', ydim='yt_ocean',
                zdim='st_ocean')
df = df.to_dataset(name='salt')
df = df.rename({'xt_ocean': 'lon', 'yt_ocean': 'lat'})
df = df.rename({'st_ocean':'depth'})
# for regridding
salt = np.copy(df.salt)
dssalt = xr.DataArray(salt, coords=[df.time, df.depth, df.lat, df.lon],
                            dims=["time","depth","lat","lon"])
dssalt = dssalt.ffill(dim='depth')
df['salt'] = dssalt
#apply regridder
dr_out = regridder(df['salt'])
varsalt = np.copy(dr_out)

time = 17.5

# Create a mosaic file
fout  = mom_dir + 'SODA_TS_IC.nc'
rg = scipy.io.netcdf_file(fout,'w')
# Dimensions
rg.createDimension('time', None)
rg.createDimension('z_l',nk)
rg.createDimension('nxp',ni)
rg.createDimension('nyp',nj)
# Variables
hnx = rg.createVariable('nxp', 'int32', ('nxp',))
hny = rg.createVariable('nyp', 'int32', ('nyp',))
hx = rg.createVariable('lon','float32',('nyp','nxp',))
hx.units = 'degrees'
hy = rg.createVariable('lat','float32',('nyp','nxp',))
hy.units = 'degrees'
hz = rg.createVariable('z_l','float32',('z_l',))
hz.units = 'meters'
tempvar  = rg.createVariable('temp','float32',('time','z_l','nyp','nxp',))
tempvar.units = 'celcius'
tempvar.missing_val = 1e20
tempvar._FillValue = 1e20
saltvar  = rg.createVariable('salt','float32',('time','z_l','nyp','nxp',))
saltvar.units = 'psu'
saltvar.missing_val = 1e20
saltvar._FillValue = 1e20
htime = rg.createVariable('time','float32',('time',))
# htime = rg.createVariable('time', 'int', ('time'))
htime.units = 'days since 1980-01-01 00:00:00'
# Values
hx[:] = lon_rho
hy[:] = lat_rho
hz[:] = np.copy(df.depth)
tempvar[:] = var
saltvar[:] = varsalt
hnx[:] = np.arange(0,ni)
hny[:] = np.arange(0,nj)
htime = time
rg.close()
