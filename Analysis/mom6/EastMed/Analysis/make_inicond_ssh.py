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
import scipy.interpolate as interp

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

# initial condition
ssh = np.zeros((1,nj,ni))

# possible tsunami
ssh[0,638,1134] = 4
ssh[0,637,1135] = 4
ssh[0,636,1136] = 4
ssh[0,635,1137] = 4
ssh[0,634,1138] = 4
ssh[0,637,1134] = 3
ssh[0,636,1135] = 3
ssh[0,635,1136] = 3
ssh[0,634,1137] = 3
ssh[0,636,1134] = 2
ssh[0,635,1135] = 2
ssh[0,634,1136] = 2
ssh[0,635,1134] = 1
ssh[0,634,1135] = 1

ssh[0,639,1135] = -4
ssh[0,638,1136] = -4
ssh[0,637,1137] = -4
ssh[0,636,1138] = -4
ssh[0,635,1139] = -4
ssh[0,639,1136] = -3
ssh[0,638,1137] = -3
ssh[0,637,1138] = -3
ssh[0,636,1139] = -3
ssh[0,639,1137] = -2
ssh[0,638,1138] = -2
ssh[0,637,1139] = -2
ssh[0,639,1138] = -1
ssh[0,638,1139] = -1

time = 17.5

# Create a mosaic file
fout = mom_dir + 'ssh_IC_EastMed.nc'
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

