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
prenameu = '3DArcticOcean_monthly_UTHMASS_'
prenamev = '3DArcticOcean_monthly_VTHMASS_'

fname = datadir + '/' + prenameu +'*.nc'
list=sorted(glob.glob(fname))
# df = xr.open_mfdataset(list)
dfu = xr.open_mfdataset(list,combine='by_coords')

fname = datadir + '/' + prenamev +'*.nc'
list=sorted(glob.glob(fname))
dfv = xr.open_mfdataset(list,combine='by_coords')

arr = np.array([datetime.datetime(1992, 1, 1) + datetime.timedelta(days=i-1) for i in range(30,365*25+7,30)])
arr = arr[:300]
newtime = xr.Dataset({'time': arr})

aa = np.arange(384,526)
bb = np.arange(386,386+71+1)

y = xr.DataArray(bb[:-1], dims='points')
x = xr.DataArray(aa[::2], dims='points')
x2 = xr.DataArray(aa[1::2], dims='points')

ds0new = xr.merge([dfv, gr])
ds = ds0new.VTHMASS*ds0new.dxG*ds0new.drF
ds = ds.to_dataset(name='transport')
vt1 = ds.transport.isel(i=x,j_g=y)
vt2 = ds.transport.isel(i=x2,j_g=y)

ds1new = xr.merge([dfu, gr])
ds1 = ds1new.UTHMASS*ds1new.dyG*ds1new.drF
ds1 = ds1.to_dataset(name='transport')
ut1 = ds1.transport.isel(i_g=x,j=y)
ut2 = ds1.transport.isel(i_g=x2,j=y)

vt1 = vt1.to_dataset(name='v1_heat_transport')
vt2 = vt2.to_dataset(name='v2_heat_transport')
ut1 = ut1.to_dataset(name='u1_heat_transport')
ut2 = ut2.to_dataset(name='u2_heat_transport')

VT_barents1 = vt1.sum(axis=(1,2))
VT_barents2 = vt2.sum(axis=(1,2))
UT_barents1 = ut1.sum(axis=(1,2))
UT_barents2 = ut2.sum(axis=(1,2))


VT_barents1.to_netcdf('dnm1.nc')
VT_barents2.to_netcdf('dnm2.nc')
UT_barents1.to_netcdf('dnm3.nc')
UT_barents2.to_netcdf('dnm4.nc')


cmd1 = 'ncks -A dnm1.nc dnm2.nc'
os.system(cmd1)
cmd2 = 'ncks -A dnm2.nc dnm3.nc'
os.system(cmd2)
cmd3 = 'ncks -A dnm3.nc dnm4.nc'
os.system(cmd3)

fname = root_folder + 'ncfiles/' + expid + '_Barents_heat_transport.nc'
cmd4 = 'cp dnm4.nc ' + fname
os.system(cmd4)

