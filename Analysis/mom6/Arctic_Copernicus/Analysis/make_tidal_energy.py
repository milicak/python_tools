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

variable = 'tideamp'
variable2 = 'h2'

mask = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/ocean_mask.nc')
df = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/OM4_025/mosaic.v20170622.unpacked/ocean_hgrid.nc')
lon_rho = np.copy(df['x'][1::2,1::2])
lat_rho = np.copy(df['y'][1::2,1::2])
nj,ni = lon_rho.shape
ds = df['x'][1::2,1::2]
ds = ds.to_dataset(name='lon')
ds['lat']=df['y'][1::2,1::2]
ds = ds.rename_dims({'nxp': 'x','nyp': 'y'})
dstmp = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/OM4_025/INPUT/tidal_amplitude.v20140616.nc')
dstmp = dstmp.rename({'ny': 'y', 'nx': 'x'})
ds[variable] = dstmp[variable]
dstmp2 = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/OM4_025/mosaic.v20170622.unpacked/ocean_topog.nc')
dstmp2 = dstmp2.rename({'ny': 'y', 'nx': 'x'})
ds[variable2] = dstmp2[variable2]

df2 = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/ocean_hgrid.nc')
lon_rho = np.copy(df2['x'][1::2,1::2])
lat_rho = np.copy(df2['y'][1::2,1::2])
nj,ni = lon_rho.shape
ds2 = df2['x'][1::2,1::2]
ds2 = ds2.to_dataset(name='lon')
ds2['lat']=df2['y'][1::2,1::2]
ds2 = ds2.rename_dims({'nxp': 'x','nyp': 'y'})

# build regridder
# regridder = xe.Regridder(ds, ds2, 'nearest_s2d')
# regridder = xe.Regridder(ds, ds2, 'nearest_s2d', reuse_weights=True)
regridder = xe.Regridder(ds, ds2, 'bilinear')
regridder = xe.Regridder(ds, ds2, 'bilinear', reuse_weights=True)

#apply regridder
dr_out = regridder(ds[variable])
dr_out = dr_out.fillna(0)
tideamp = np.copy(dr_out)
dr_out2 = regridder(ds[variable2])
dr_out2 = dr_out2.fillna(0)
h2 = np.copy(dr_out2)
nx = np.copy(lon_rho[0,:])
ny = np.copy(lat_rho[:,0])

# Create a mosaic file
rg = scipy.io.netcdf_file('tidal_amplitude_Arctic.nc','w')
# Dimensions
rg.createDimension('nx',ni)
rg.createDimension('ny',nj)
# Variables
hx = rg.createVariable('nx','float32',('nx',))
hx.units = 'degrees east'
hx.cartesian_axis = 'X'
hy = rg.createVariable('ny','float32',('ny',))
hy.units = 'degrees north'
hy.cartesian_axis = 'Y'
tidevar  = rg.createVariable('tideamp','float32',('ny','nx',))
tidevar.units = 'm s-1'
tidevar.missing_val = 1e20
tidevar._FillValue = 1e20
h2var  = rg.createVariable('h2','float32',('ny','nx',))
h2var.units = 'meters^2'
h2var.standard_name = 'Variance of sub-grid scale topography'
# Values
hx[:] = nx
hy[:] = ny
h2var[:] = h2
tidevar[:] = tideamp
rg.close()



