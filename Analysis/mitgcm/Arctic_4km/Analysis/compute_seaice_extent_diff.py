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

prename = '2DArcticOcean_'

datadir0 = root_folder+expid0
fname = datadir0+'/'+prename+'*'+'SIarea*'
list=sorted(glob.glob(fname))
df0 = xr.open_mfdataset(list)

datadir1 = root_folder+expid1
fname = datadir1+'/'+prename+'*'+'SIarea*'
list=sorted(glob.glob(fname))
df1 = xr.open_mfdataset(list)

datadir3 = root_folder+expid3
fname = datadir3+'/'+prename+'*'+'SIarea*'
list=sorted(glob.glob(fname))
df3 = xr.open_mfdataset(list)

# time = pd.date_range('1992-01-01', freq='M', periods=12 * 25)
# df['time'] = time

# difference in sea ice extent
df = df1.SIarea-df0.SIarea
# Compute anomalies by removing the time-mean.
df = df - df.mean(dim='time')
df = df.to_dataset(name='SIarea')
df.to_netcdf('/archive/milicak/MITgcm_c65/Projects/Arctic_4km/ncfiles/Exp02_1_seaice_area_diff.nc')

df = df3.SIarea-df0.SIarea
# Compute anomalies by removing the time-mean.
df = df - df.mean(dim='time')
df = df.to_dataset(name='SIarea')
df.to_netcdf('/archive/milicak/MITgcm_c65/Projects/Arctic_4km/ncfiles/Exp02_3_seaice_area_diff.nc')

