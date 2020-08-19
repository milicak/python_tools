import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import os
from datetime import date
import datetime

root_folder = '/archive/milicak/MITgcm_c65/Projects/Arctic_4km/'

# expid = 'Exp02_0';
# expid = 'Exp02_1';
# expid = 'Exp02_2';
expid = 'Exp02_3';
print(expid)

gridname = root_folder + expid + '/grid.nc'
gr = xr.open_dataset(gridname)
gr = gr.drop_dims({'time'})

fyear = 1992
# fyear = 1992
lyear = 2018
# lyear = fyear+1
datadir = root_folder+expid
# os.chdir(datadir)
prenameu = '3DArcticOcean_monthly_UVELMASS_'

fname = datadir + '/' + prenameu +'*.nc'
list=sorted(glob.glob(fname))
dfu = xr.open_mfdataset(list,combine='by_coords')

aa = np.repeat(np.array([1449]),30)
bb = np.arange(659,689)

x = xr.DataArray(aa, dims='points')
y = xr.DataArray(bb, dims='points')

ds1new = xr.merge([dfu, gr])
ds1 = ds1new.UVELMASS*ds1new.dyG*ds1new.drF

# vt = ds1.transport[:,:,659:688,1449]
# VT_bering = vt.sum(axis=(1,2))
# dnm = VT_bering.compute()

ds1 = ds1.to_dataset(name='transport')
ut1 = ds1.transport.isel(i_g=x,j=y)

ut1 = ut1.to_dataset(name='volume_transport')
UT_bering = ut1.sum(axis=(1,2))

fname = root_folder + 'ncfiles/' + expid + '_Bering_volume_transport.nc'
UT_bering.to_netcdf(fname)






