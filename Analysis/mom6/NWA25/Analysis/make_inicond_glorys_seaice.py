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

variable = 'aice'
variable1 = 'hi'
variable2 = 'siconc'
variable3 = 'sithick'

root_folder = '/Volumes/A1/workdir/milicak/datasets/MOM6/NWA25/'
mom_dir = '/Volumes/A1/workdir/milicak/datasets/MOM6/NWA25/'
path_regional_grid = mom_dir + './ocean_hgrid.nc'                         

df2 = xr.open_dataset(path_regional_grid)                                       
lon_rho = np.copy(df2['x'][1::2,1::2])                                          
lat_rho = np.copy(df2['y'][1::2,1::2])                                          
nj,ni = lon_rho.shape                                                           
ds2 = xr.Dataset()                                                              
ds2['lon'] = df2['x'][1::2,1::2]                                                
ds2['lat'] = df2['y'][1::2,1::2]                                                
                                                                                
# ice concentration should be between 0 and 1
df = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/glorys/tmp/GLORYS_REANALYSIS_1996-01-01.nc')
df2['lon'] = df['longitude']                                                    
df2['lat'] = df['latitude']                                                     

# build regridder
# regridder = xe.Regridder(df, ds2, 'nearest_s2d')
regridder = xe.Regridder(df2, ds2, 'bilinear', reuse_weights=True)

#apply regridder
dr_out = regridder(df[variable2])
dr_out = dr_out.fillna(0)
dr_out2 = regridder(df[variable3])
dr_out2 = dr_out2.fillna(0)
aice = np.copy(dr_out)
hice = np.copy(dr_out2)

aice[aice>1] = 1

time = 17.5

# Create a mosaic file
fout = mom_dir + 'glorys_ic_seaice_1996.nc'
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
aicevar  = rg.createVariable('aice','float32',('time','nyp','nxp',))
aicevar.units = 'concentration 0-1'
aicevar.missing_val = 1e20
aicevar._FillValue = 1e20
hicevar  = rg.createVariable('hice','float32',('time','nyp','nxp',))
hicevar.units = 'meters'
hicevar.missing_val = 1e20
hicevar._FillValue = 1e20
htime = rg.createVariable('time','float32',('time',))
# htime = rg.createVariable('time', 'int', ('time'))
htime.units = 'days since 1980-01-01 00:00:00'
# Values
hx[:] = lon_rho
hy[:] = lat_rho
aicevar[:] = aice
hicevar[:] = hice
htime = time
rg.close()

