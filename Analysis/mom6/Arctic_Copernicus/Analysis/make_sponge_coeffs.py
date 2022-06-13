import os
import numpy as np
import xesmf as xe
import xarray as xr
import scipy.io
from scipy.io import savemat
from scipy.io import loadmat
from scipy.signal import medfilt2d
import netCDF4
from scipy.interpolate import griddata


df2 = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/ocean_hgrid.nc')
lon_rho = np.copy(df2['x'][1::2,1::2])
lat_rho = np.copy(df2['y'][1::2,1::2])
nj,ni = lon_rho.shape
X = np.int32(np.arange(1,ni+1))
Y = np.int32(np.arange(1,nj+1))

# close to boundary time scale
time1  = 2.5 # 5 days
# away from boundary time scale
time2  = 20 # 30 days
# total grid points to sponge
grdpoints = 50

coeffs=1/np.linspace(time1,time2,grdpoints)/86400

sponge_var1 = np.zeros((nj,ni))
tmp = np.tile(coeffs,(ni,1))
sponge_var1[1:grdpoints+1,:] = np.transpose(tmp)
sponge_var1[-grdpoints-1:-1,:] = np.flip(np.transpose(tmp))

sponge_var2 = np.zeros((nj,ni))
tmp = np.tile(coeffs,(nj,1))
sponge_var2[:,1:grdpoints+1] = tmp
sponge_var2[:,-grdpoints-1:-1] = np.flip(tmp)

sponge_var = np.zeros((nj,ni))
sponge_var = np.maximum(sponge_var1,sponge_var2)
sponge_var[0,:] = 0
sponge_var[-1,:] = 0
sponge_var[:,0] = 0
sponge_var[:,-1] = 0

# Create a mosaic file
rg = scipy.io.netcdf_file('ocean_sponge.nc','w')
# Dimensions
rg.createDimension('X',ni)
rg.createDimension('Y',nj)
# Variables
hx = rg.createVariable('X','i4',('X',))
hx.point_spacing = 'even'
hx.axis = 'X'
hy = rg.createVariable('Y','i4',('Y',))
hy.point_spacing = 'even'
hy.axis = 'Y'
sponge  = rg.createVariable('IDAMP','float32',('Y','X',))
sponge.long_name = 'Damping coefficient'
sponge.units = 's-1'
sponge._FillValue = 1e20
sponge[:] = sponge_var
hx[:] = X
hy[:] = Y
rg.close()

