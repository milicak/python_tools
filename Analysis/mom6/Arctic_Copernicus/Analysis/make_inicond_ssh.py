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

variable = 'sealv'

ds = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic/NOIIAJRAOC20TR_TL319_tn14_20190710.micom.hm.1958-01.nc')
gr1 = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic/NorESM_grid.nc')
mask = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/ocean_mask.nc')

df = flood_kara(ds[variable], xdim='x', ydim='y', zdim='z')
df = df[0,0,:,:]
df = df.to_dataset(name='zeta')
df['lon'] = gr1['plon']
df['lat'] = gr1['plat']

df2 = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/ocean_hgrid.nc')
lon_rho = np.copy(df2['x'][1::2,1::2])
lat_rho = np.copy(df2['y'][1::2,1::2])
nj,ni = lon_rho.shape
ds2 = df2['x'][1::2,1::2]
ds2 = ds2.to_dataset(name='lon')
ds2['lat']=df2['y'][1::2,1::2]
ds2 = ds2.rename_dims({'nxp': 'x','nyp': 'y'})

zeta = np.copy(df.zeta)
dszeta = xr.DataArray(zeta, coords=[df.y, df.x],
                            dims=["y","x"])
df['zeta'] = dszeta

# build regridder
regridder = xe.Regridder(df, ds2, 'nearest_s2d')
regridder = xe.Regridder(df, ds2, 'nearest_s2d', reuse_weights=True)

#apply regridder
dr_out = regridder(df['zeta'])
zeta = np.copy(dr_out)*np.copy(mask.mask)

# create 3d ssh variable
# not sure if i need to multiply with mask
ssh = np.reshape(zeta,(1,zeta.shape[0],zeta.shape[1]))


time = 17.5

# Create a mosaic file
rg = scipy.io.netcdf_file('ssh_IC.nc','w')
# Dimensions
rg.createDimension('time', None)
rg.createDimension('nxp',ni)
rg.createDimension('nyp',nj)
# Variables
hx = rg.createVariable('lon','float32',('nyp','nxp',))
hx.units = 'degrees'
hy = rg.createVariable('lat','float32',('nyp','nxp',))
hy.units = 'degrees'
sshvar  = rg.createVariable('ssh','float32',('time','nyp','nxp',))
sshvar.units = 'meters'
sshvar.missing_val = 1e20
sshvar._FillValue = 1e20
htime = rg.createVariable('time','float32',('time',))
# htime = rg.createVariable('time', 'int', ('time'))
htime.units = 'days since 1958-01-01 00:00:00'
# Values
hx[:] = lon_rho
hy[:] = lat_rho
sshvar[:] = ssh
htime = time
rg.close()



