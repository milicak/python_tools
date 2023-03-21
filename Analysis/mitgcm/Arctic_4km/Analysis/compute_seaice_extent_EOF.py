import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import os
from datetime import date
from eofs.xarray import Eof

root_folder = '/archive/milicak/MITgcm_c65/Projects/Arctic_4km/'

expid0 = 'Exp02_0';
# Atlantic v1
expid1 = 'Exp02_1';
# Atlantic v2
expid3 = 'Exp02_3';

gridname = root_folder + expid0 + '/grid.nc'
gr = xr.open_dataset(gridname)

df1 = xr.open_dataset('/archive/milicak/MITgcm_c65/Projects/Arctic_4km/ncfiles/Exp02_1_seaice_area_diff.nc')
df3 = xr.open_dataset('/archive/milicak/MITgcm_c65/Projects/Arctic_4km/ncfiles/Exp02_3_seaice_area_diff.nc')


# Create an EOF solver to do the EOF analysis. Square-root of cosine of
# latitude weights are applied before the computation of EOFs.
coslat = np.cos(np.deg2rad(gr['YC'].values)).clip(0., 1.)
coslat = np.copy(coslat[:,855])
wgts = np.sqrt(coslat)[..., np.newaxis]
wgts = np.ones((1536,1))
solver1 = Eof(df1.SIarea, weights=wgts)
# Retrieve the leading EOF, expressed as the covariance between the leading PC
# time series and the input SLP anomalies at each grid point.
# eof1 = solver1.eofsAsCovariance(neofs=1)
eof1 = solver1.eofsAsCorrelation(neofs=1)
ds = eof1.to_dataset(name='eof1')
ds.to_netcdf('/archive/milicak/MITgcm_c65/Projects/Arctic_4km/ncfiles/Exp02_1_seaice_area_diff_EOF.nc')
# pc1 = solver1.pcs(npcs=1, pcscaling=1)
variance_fraction1 = solver1.varianceFraction(neigs=1)

solver3 = Eof(df3.SIarea, weights=wgts)
eof3 = solver3.eofsAsCorrelation(neofs=1)
# pc3 = solver3.pcs(npcs=1, pcscaling=1)
variance_fraction3 = solver3.varianceFraction(neigs=1)
ds3 = eof3.to_dataset(name='eof1')
ds3.to_netcdf('/archive/milicak/MITgcm_c65/Projects/Arctic_4km/ncfiles/Exp02_3_seaice_area_diff_EOF.nc')


