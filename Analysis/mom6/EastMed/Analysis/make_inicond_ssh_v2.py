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

mom_dir = '/okyanus/users/milicak/dataset/MOM6/EastMed/'
path_regional_grid = mom_dir + './ocean_hgrid.nc'

df2 = xr.open_dataset(path_regional_grid)
lon_rho = np.copy(df2['x'][1::2,1::2])
lat_rho = np.copy(df2['y'][1::2,1::2])
nj,ni = lon_rho.shape
ds2 = df2['x'][1::2,1::2]
ds2 = ds2.to_dataset(name='lon')
ds2['lat']=df2['y'][1::2,1::2]
ds2 = ds2.rename_dims({'nxp': 'x','nyp': 'y'})


# x = np.linspace(1050,1150,3)
# x2 = np.linspace(1150,1250,8)
# y = np.linspace(420,465,3)
# y2 = np.linspace(465,560,8)
# x = np.concatenate((x,x2))
# y = np.concatenate((y,y2))
# x = np.int16(x)
# y = np.int16(y)

coeff = 0.25
# initial condition
ssh = np.zeros((1,nj,ni))
# possible tsunami
x = np.arange(1050,1150,2)
y = np.arange(420,470)
for itr in range(0,x.size):
    ssh[0,y[itr],x[itr]] = 4*coeff
    ssh[0,y[itr]+1,x[itr]] = 3*coeff
    ssh[0,y[itr]+2,x[itr]] = 2*coeff
    ssh[0,y[itr]+3,x[itr]] = 1*coeff
    ssh[0,y[itr]-1,x[itr]] = -4*coeff
    ssh[0,y[itr]-2,x[itr]] = -3*coeff
    ssh[0,y[itr]-3,x[itr]] = -2*coeff
    ssh[0,y[itr]-4,x[itr]] = -1*coeff

x = np.arange(1051,1150,2)
for itr in range(0,x.size):
    ssh[0,y[itr],x[itr]] = 4*coeff
    ssh[0,y[itr]+1,x[itr]] = 3*coeff
    ssh[0,y[itr]+2,x[itr]] = 2*coeff
    ssh[0,y[itr]+3,x[itr]] = 1*coeff
    ssh[0,y[itr]-1,x[itr]] = -4*coeff
    ssh[0,y[itr]-2,x[itr]] = -3*coeff
    ssh[0,y[itr]-3,x[itr]] = -2*coeff
    ssh[0,y[itr]-4,x[itr]] = -1*coeff

x = np.arange(1150,1220)
y = np.arange(470,540)
for itr in range(0,x.size):
    ssh[0,y[itr],x[itr]] = 4*coeff
    ssh[0,y[itr]+1,x[itr]] = 3*coeff
    ssh[0,y[itr]+2,x[itr]] = 2*coeff
    ssh[0,y[itr]+3,x[itr]] = 1*coeff
    ssh[0,y[itr]-1,x[itr]] = -4*coeff
    ssh[0,y[itr]-2,x[itr]] = -3*coeff
    ssh[0,y[itr]-3,x[itr]] = -2*coeff
    ssh[0,y[itr]-4,x[itr]] = -1*coeff

ssh = -ssh
time = 17.5

# Create a mosaic file
fout = mom_dir + 'ssh_IC_EastMed_v2.nc'
rg = scipy.io.netcdf_file(fout,'w')
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
htime.units = 'days since 1980-01-01 00:00:00'
# Values
hx[:] = lon_rho
hy[:] = lat_rho
sshvar[:] = ssh
htime = time
rg.close()

