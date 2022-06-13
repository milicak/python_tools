import os
import numpy as np
import xarray as xr
import scipy.io
import netCDF4
import matplotlib.pyplot as plt
import xesmf
from HCtFlood.kara import flood_kara


root_folder = '/okyanus/users/milicak/dataset/MOM6/TSS/'
fname1 = 'REG_0_1_0_salinity_0.01.nc'
ds = xr.open_dataset(root_folder + fname1)
fname1 = 'REG_0_1_0_temperature_0.01.nc'
ds1 = xr.open_dataset(root_folder + fname1)
fname1 = 'REG_0_1_0_water_level_0.01.nc'
ds2 = xr.open_dataset(root_folder + fname1)
ds = xr.merge([ds,ds1,ds2])
ds = ds.where(ds>=4)

mom_dir = '/okyanus/users/milicak/dataset/MOM6/TSS/'
path_regional_grid = mom_dir + './ocean_hgrid.nc'

ds = ds.where(ds>=4)
ds = ds.ffill('level')

ds['temperature'][0,-1,:,:] = ds['temperature'][0,-2,:,:]
ds['salinity'][0,-1,:,:] = ds['salinity'][0,-2,:,:]

dft = flood_kara(ds['temperature'], xdim='lon', ydim='lat', zdim='level')
dst = flood_kara(ds['salinity'], xdim='lon', ydim='lat', zdim='level')
dzt = flood_kara(ds['water_level'], xdim='lon', ydim='lat')

vartemp = np.copy(dft)
nk = vartemp.shape[1]
nj = vartemp.shape[2]
ni = vartemp.shape[3]
varsalt = np.copy(dst)
varsalt += 0.15
# varsalt[varsalt==0.15] = 0.0
ssh = np.copy(dzt[:,0,:,:])

time = 17.5

# Create a mosaic file
fout  = mom_dir + 'shyfem_TSssh_ICnew.nc'
rg = scipy.io.netcdf_file(fout,'w')
# Dimensions
rg.createDimension('time', None)
rg.createDimension('depth',nk)
rg.createDimension('longitude',ni)
rg.createDimension('latitude',nj)
# Variables
hnx = rg.createVariable('longitude', 'float32', ('longitude',))
hnx.units = 'degrees east'
hnx.standard_name = 'longitude'
hny = rg.createVariable('latitude', 'float32', ('latitude',))
hny.units = 'degrees north'
hny.standard_name = 'latitude'
hz = rg.createVariable('depth','float32',('depth',))
hz.units = 'meters'
hz._CoordinateZisPositive = 'down'
tempvar  = rg.createVariable('temp','float32',('time','depth','latitude','longitude',))
tempvar.units = 'celcius'
tempvar.missing_val = 1e20
tempvar._FillValue = 1e20
saltvar  = rg.createVariable('salt','float32',('time','depth','latitude','longitude',))
saltvar.units = 'psu'
saltvar.missing_val = 1e20
saltvar._FillValue = 1e20
sshvar  = rg.createVariable('ssh','float32',('time','latitude','longitude',))
sshvar.units = 'meters'
sshvar.missing_val = 1e20
sshvar._FillValue = 1e20
htime = rg.createVariable('time','float32',('time',))
# htime = rg.createVariable('time', 'int', ('time'))
htime.units = 'days since 1996-01-01 00:00:00'
# Values
hnx[:] = np.copy(ds.lon)
hny[:] = np.copy(ds.lat)
hz[:] = np.copy(ds.level)
tempvar[:] = vartemp
saltvar[:] = varsalt
sshvar[:] = ssh
htime = time
rg.close()

# Create a mosaic file
fout  = mom_dir + 'vgrid_93_1m.nc'
rg = scipy.io.netcdf_file(fout,'w')
dz = np.concatenate(([1],ds.level[1:].data-ds.level[:-1].data))
rg.createDimension('nz',nk)
hz = rg.createVariable('dz','double',('nz',))
hz.units = 'm'
hz.long_name = 'z coordinate level thickness'
hz[:] = np.copy(dz)
rg.close()

