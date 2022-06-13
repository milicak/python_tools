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

df2 = xr.open_dataset(path_regional_grid)
lon_rho = np.copy(df2['x'][1::2,1::2])
lat_rho = np.copy(df2['y'][1::2,1::2])
nj,ni = lon_rho.shape

target_grid = xr.open_dataset(path_regional_grid)
# target_grid['x'] -= 360.0
target_t = (
   target_grid
   [['x', 'y']]
   .isel(nxp1=slice(1, None, 2), nyp1=slice(1, None, 2))
   .rename({'y': 'lat', 'x': 'lon', 'nxp1': 'xh', 'nyp1': 'yh'})
)


regrid_kws = dict(method='nearest_s2d', reuse_weights=False, periodic=False)
regrid_kws = dict(method='nearest_s2d', reuse_weights=True, periodic=False)
shyfem_to_t = xesmf.Regridder(ds, target_t, filename='regrid_shyfem_tracers.nc', **regrid_kws)
df = shyfem_to_t(ds)

dft = flood_kara(df['temperature'], xdim='xh', ydim='yh', zdim='level')
dst = flood_kara(df['salinity'], xdim='xh', ydim='yh', zdim='level')
dzt = flood_kara(df['water_level'], xdim='xh', ydim='yh')

vartemp = np.copy(dft)
nk = vartemp.shape[1]
varsalt = np.copy(dst)
varsalt += 0.15
varsalt[varsalt==0.15] = 0.0
ssh = np.copy(dzt[:,0,:,:])

time = 17.5

# Create a mosaic file
fout  = mom_dir + 'shyfem_TSssh_IC.nc'
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
hz[:] = np.copy(ds.level)
tempvar[:] = vartemp
saltvar[:] = varsalt
sshvar[:] = ssh
hnx[:] = np.arange(0,ni)
hny[:] = np.arange(0,nj)
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

