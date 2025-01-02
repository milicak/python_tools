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

variable = 'sos'

fname1 = '/okyanus/users/milicak/dataset/MOM6/OM4_025/INPUT.JRA.v2019.07.04/salt_restore_JRA.1440x1080.v20190706.nc'
ds = xr.open_dataset(fname1)
mom_dir = '~/dataset/MOM6/Arctic_Copernicus/'
path_regional_grid = mom_dir + './ocean_hgrid.nc'

df2 = xr.open_dataset(path_regional_grid)
lon_rho = np.copy(df2['x'][1::2,1::2])
lat_rho = np.copy(df2['y'][1::2,1::2])
nj,ni = lon_rho.shape
ds2 = df2['x'][1::2,1::2]
ds2 = ds2.to_dataset(name='lon')
ds2['lat']=df2['y'][1::2,1::2]
ds2 = ds2.rename_dims({'nxp': 'x','nyp': 'y'})

dft = flood_kara(ds['sos'], xdim='i', ydim='j')
dft.load();

# build regridder
regridder = xe.Regridder(ds, ds2, 'patch', reuse_weights=True)
regridder = xe.Regridder(ds, ds2, 'patch', reuse_weights=True, periodic=True)

#apply regridder
dr_out = regridder(dft[:,0,:,:])
dr_out1 = dr_out.where(dr_out!=0)
dr_out1 = dr_out1.interpolate_na(dim="y", method="nearest",
                         fill_value="extrapolate")
dfs = dr_out1.to_dataset(name='sss')

# Create a mosaic file
fout = mom_dir + 'SSS_climatology.nc'
dfs.to_netcdf(fout,unlimited_dims='time')
