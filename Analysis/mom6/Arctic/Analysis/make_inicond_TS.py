import os
import numpy as np
import xesmf as xe
import xarray as xr
import scipy.io
from scipy.io import savemat
from scipy.io import loadmat
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.colors as colors
from scipy.signal import medfilt2d
import netCDF4
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.path import Path
#for interpolation
from scipy.spatial import cKDTree
from HCtFlood.kara import flood_kara

variable = 'templvl'

ds = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic/NOIIAJRAOC20TR_TL319_tn14_20190710.micom.hm.1958-01.nc')
gr1 = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic/NorESM_grid.nc')
gr2 = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic/ocean_geometry.nc')
mask = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic/ocean_mask.nc')

df = flood_kara(ds[variable], xdim='x', ydim='y', zdim='depth')
# df = df[0,0,:,:]
df = df.to_dataset(name='temp')
df['lon'] = gr1['plon']
df['lat'] = gr1['plat']

# output grid info
df2 = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic/ocean_hgrid.nc')
lon_rho = np.copy(df2['x'][1::2,1::2])
lat_rho = np.copy(df2['y'][1::2,1::2])
nj,ni = lon_rho.shape
ds2 = df2['x'][1::2,1::2]
ds2 = ds2.to_dataset(name='lon')
ds2['lat']=df2['y'][1::2,1::2]
ds2 = ds2.rename_dims({'nxp1': 'x','nyp1': 'y'})

# for regridding
temp = np.copy(df.temp)
dstemp = xr.DataArray(temp, coords=[df.time, df.depth, df.y, df.x],
                            dims=["time","depth","y","x"])
df['temp'] = dstemp

# build regridder
regridder = xe.Regridder(df, ds2, 'nearest_s2d', reuse_weights=True)

#apply regridder
dr_out = regridder(df['temp'])

# create 3d ssh variable
# not sure if i need to multiply with mask
var = np.copy(dr_out)
nk = var.shape[1]


# salinity
variable = 'salnlvl'
df = flood_kara(ds[variable], xdim='x', ydim='y', zdim='depth')
df = df.to_dataset(name='salt')
df['lon'] = gr1['plon']
df['lat'] = gr1['plat']
# for regridding
salt = np.copy(df.salt)
dssalt = xr.DataArray(salt, coords=[df.time, df.depth, df.y, df.x],
                            dims=["time","depth","y","x"])
df['salt'] = dssalt
#apply regridder
dr_out = regridder(df['salt'])
varsalt = np.copy(dr_out)


time = 17.5

# Create a mosaic file
rg = scipy.io.netcdf_file('tempsalt_IC.nc','w')
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
htime.units = 'days since 1958-01-01 00:00:00'
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



