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
import xesmf as xe

variable = 'chlor_a'

fname1 = '/Volumes/A1/workdir/milicak/datasets/MOM6/CORE2/NYF_v2.0/seawifs-clim-1997-2010.nc'
ds = xr.open_dataset(fname1)
mom_dir = '/Volumes/A1/workdir/milicak/datasets/MOM6/NWA25_v2/'
path_regional_grid = mom_dir + './ocean_hgrid.nc'

df2 = xr.open_dataset(path_regional_grid)
lon_rho = np.copy(df2['x'][1::2,1::2])
lat_rho = np.copy(df2['y'][1::2,1::2])
nj,ni = lon_rho.shape
ds2 = df2['x'][1::2,1::2]
ds2 = ds2.to_dataset(name='lon')
ds2['lat']=df2['y'][1::2,1::2]
ds2 = ds2.rename_dims({'nxp': 'x','nyp': 'y'})

dft = flood_kara(ds['chlor_a'], xdim='lon', ydim='lat')

# build regridder
regridder = xe.Regridder(ds, ds2, 'bilinear', reuse_weights=True)

#apply regridder
dr_out = regridder(dft[:,0,:,:])
dfs = dr_out.to_dataset(name=variable)

# Create a mosaic file
fout = mom_dir + 'chl_climatology.nc'
dfs.to_netcdf(fout,unlimited_dims='time')  
