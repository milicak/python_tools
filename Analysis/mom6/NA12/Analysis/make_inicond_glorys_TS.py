import os
import numpy as np
import xarray as xr
import scipy.io
import netCDF4
import matplotlib.pyplot as plt


root_folder = '/okyanus/users/milicak/dataset/MOM6/NA12/'
fname1 = 'out2.nc'
ds = xr.open_dataset(root_folder + fname1)
mom_dir = '/okyanus/users/milicak/dataset/MOM6/NA12/'
path_regional_grid = mom_dir + './ocean_hgrid.nc'

df2 = xr.open_dataset(path_regional_grid)
lon_rho = np.copy(df2['x'][1::2,1::2])
lat_rho = np.copy(df2['y'][1::2,1::2])
nj,ni = lon_rho.shape


vartemp = np.copy(ds.temp)
nk = vartemp.shape[1]
varsalt = np.copy(ds.salt)
ssh = np.copy(ds.ssh)

time = 17.5

# Create a mosaic file
fout  = mom_dir + 'glorys_TS_IC.nc'
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
sshvar  = rg.createVariable('ssh','float32',('time','nyp','nxp',))
sshvar.units = 'meters'
sshvar.missing_val = 1e20
sshvar._FillValue = 1e20
htime = rg.createVariable('time','float32',('time',))
# htime = rg.createVariable('time', 'int', ('time'))
htime.units = 'days since 1996-01-01 00:00:00'
# Values
hx[:] = lon_rho
hy[:] = lat_rho
hz[:] = np.copy(ds.z_l)
tempvar[:] = vartemp
saltvar[:] = varsalt
sshvar[:] = ssh
hnx[:] = np.arange(0,ni)
hny[:] = np.arange(0,nj)
htime = time
rg.close()

