# This file is designed to be cut and pasted into an ipython --pylab
# session. Otherwise, you'll need to "import np as np" then
# convert "array" to "np.array".
import os
import numpy as np
import xarray as xr
import scipy.io
import netCDF4
import xesmf


ds = xr.open_dataset('~/dataset/world_bathy/bathy_marmara.nc')
ds = ds.rename({'y':'lat','x': 'lon'})
df = xr.open_dataset('ocean_hgrid.nc')
lon_rho = np.copy(df['x'][1::2,1::2])
lat_rho = np.copy(df['y'][1::2,1::2])

target_grid = xr.open_dataset('ocean_hgrid.nc')
target_t = (
   target_grid
   [['x', 'y']]
   .isel(nxp=slice(1, None, 2), nyp=slice(1, None, 2))
   .rename({'y': 'lat', 'x': 'lon', 'nxp': 'xh', 'nyp': 'yh'})
)

# regrid_kws = dict(method='bilinear', reuse_weights=False, periodic=False,
regrid_kws = dict(method='bilinear', reuse_weights=True, periodic=False,
                  ignore_degenerate=True)
TSS = xesmf.Regridder(ds, target_t, filename='regrid_TSS.nc', **regrid_kws)
tmp = TSS(ds.z)

