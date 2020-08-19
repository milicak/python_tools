import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import os
from datetime import date

root_folder = '/archive/milicak/MITgcm_c65/Projects/Arctic_4km/'

expid = 'Exp02_3';

gridname = root_folder + expid + '/grid.nc'
gr = xr.open_dataset(gridname)

fyear = 1992
# fyear = 1992
lyear = 2018
# lyear = fyear+1
datadir = root_folder+expid
# os.chdir(datadir)
prename = '2DArcticOcean_'
postname = '_avg.nc'

fname = datadir+'/'+prename+'*'+'SIarea*'
list=sorted(glob.glob(fname))

time = pd.date_range('1992-01-01', freq='M', periods=12 * 25)

df = xr.open_mfdataset(list)
# df['time'] = time

si = df*gr.rA
SI = si.sum(dim=['i','j'])
fname = root_folder + 'ncfiles/' + expid + '_seaice_area.nc'
SI.to_netcdf(fname)
siann = SI.groupby('time.year').mean('time')*1e-12

aa = xr.where(df.SIarea<0.15,0,1)
aa = aa*gr.rA
SIext = aa.sum(dim=['i','j'])
ds = SIext.to_dataset(name='SI_extent')
fname = root_folder + 'ncfiles/' + expid + '_seaice_extent.nc'
ds.to_netcdf(fname)

