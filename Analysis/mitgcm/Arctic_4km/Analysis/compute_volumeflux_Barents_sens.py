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
prenamev = '3DArcticOcean_monthly_VVELMASS_'

fname = datadir + '/' + prenameu +'*.nc'
list=sorted(glob.glob(fname))
# df = xr.open_mfdataset(list)
dfu = xr.open_mfdataset(list,combine='by_coords')

fname = datadir + '/' + prenamev +'*.nc'
list=sorted(glob.glob(fname))
dfv = xr.open_mfdataset(list,combine='by_coords')

# old way
# time = pd.date_range('1992-01-01', freq='M', periods=12 * 25)
# new way
# arr = np.array([datetime.datetime(1991, 12, 16) + datetime.timedelta(days=i-1) for i in range(30,365*25+7,30)])
arr = np.array([datetime.datetime(1992, 1, 1) + datetime.timedelta(days=i-1) for i in range(30,365*25+7,30)])
arr = arr[:300]
newtime = xr.Dataset({'time': arr})
# df['time'] = newtime['time']

aa = np.arange(384,526)
bb = np.arange(386,386+71+1)
# just v velocity
# aa[1::2],bb[:-1]
# both u and v velocity
# aa[::2],bb[:-1]

y = xr.DataArray(bb[:-1], dims='points')
x = xr.DataArray(aa[::2], dims='points')
x2 = xr.DataArray(aa[1::2], dims='points')

# volume transport at fram strait m^3/s
# dynew = df.dyC.rename({'j_g': 'j', 'i': 'i_g'})
ds0new = xr.merge([dfv, gr])
# ds = dfv.VVELMASS*gr.dxG*gr.drF
ds = ds0new.VVELMASS*ds0new.dxG*ds0new.drF
ds = ds.to_dataset(name='transport')
vt1 = ds.transport.isel(i=x,j_g=y)
vt2 = ds.transport.isel(i=x2,j_g=y)

# i = 0
# vt1 = ds.transport.isel(j_g=bb[i],i=aa[2*i+1])
# vt2 = ds.transport.isel(j_g=bb[i],i=aa[2*i])
# for i in range(1,bb.shape[0]-1):
#     print(i,bb[i])
#     tmp1 = ds.transport.isel(j_g=bb[i],i=aa[2*i+1])
#     tmp2 = ds.transport.isel(j_g=bb[i],i=aa[2*i])
#     vt1 = vt1 + tmp1
#     vt2 = vt2 + tmp2
#

# vt1 = np.zeros((ds.transport.shape[0],ds.transport.shape[1],bb.shape[0]-1))
# vt2 = np.zeros((ds.transport.shape[0],ds.transport.shape[1],bb.shape[0]-1))
# for i in range(0,bb.shape[0]-1):
#     print(bb[i])
#     vt1[:,:,i] = np.copy(ds.transport[:,:,bb[i],aa[2*i+1]])
#     vt2[:,:,i] = np.copy(ds.transport[:,:,bb[i],aa[2*i]])


# vt = ds.data[:,:,456,205:531]

ds1new = xr.merge([dfu, gr])
# ds1 = dfu.UVELMASS*gr.dyG*gr.drF
ds1 = ds1new.UVELMASS*ds1new.dyG*ds1new.drF
ds1 = ds1.to_dataset(name='transport')
ut1 = ds1.transport.isel(i_g=x,j=y)
ut2 = ds1.transport.isel(i_g=x2,j=y)

# i = 0
# ut1 = ds1.transport.isel(j=bb[i],i_g=aa[2*i+1])
# ut2 = ds1.transport.isel(j=bb[i],i_g=aa[2*i])
# for i in range(1,bb.shape[0]-1):
#     print(i,bb[i])
#     tmp1 = ds.transport.isel(j=bb[i],i_g=aa[2*i+1])
#     tmp2 = ds.transport.isel(j=bb[i],i_g=aa[2*i])
#     ut1 = ut1 + tmp1
#     ut2 = ut2 + tmp2
#
# ut1 = np.zeros((ds.transport.shape[0],ds.transport.shape[1],bb.shape[0]-1))
# ut2 = np.zeros((ds.transport.shape[0],ds.transport.shape[1],bb.shape[0]-1))
# for i in range(0,bb.shape[0]-1):
#     print(bb[i])
#     ut1[:,:,i] = np.copy(ds.transport[:,:,bb[i],aa[2*i+1]])
#     ut2[:,:,i] = np.copy(ds.transport[:,:,bb[i],aa[2*i]])

vt1 = vt1.to_dataset(name='v1_volume_transport')
vt2 = vt2.to_dataset(name='v2_volume_transport')
ut1 = ut1.to_dataset(name='u1_volume_transport')
ut2 = ut2.to_dataset(name='u2_volume_transport')

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

# dsnew = xr.merge([VT_barents1, VT_barents2, UT_barents1, UT_barents2])

# TT_barents = VT_barents1 + VT_barents2 + UT_barents
# dnm = TT_barents.compute()

# data = xr.DataArray(dnm, dims=('time'), coords={'time': time})
# data1 = xr.DataArray(VT_barents1, dims=('time'), coords={'time': dfv['time']})
# data1 = data1.to_dataset(name='v1_volume_transport')
# data2 = xr.DataArray(VT_barents2, dims=('time'), coords={'time': dfv['time']})
# data2 = data2.to_dataset(name='v2_volume_transport')
# data3 = xr.DataArray(UT_barents1, dims=('time'), coords={'time': dfv['time']})
# data3 = data3.to_dataset(name='u1_volume_transport')
# data4 = xr.DataArray(UT_barents2, dims=('time'), coords={'time': dfv['time']})
# data4 = data4.to_dataset(name='u2_volume_transport')
# dsnew = xr.merge([data1, data2, data3, data4])


fname = root_folder + 'ncfiles/' + expid + '_Barents_volume_transport.nc'
cmd4 = 'cp dnm4.nc ' + fname
os.system(cmd4)
# dsnew.to_netcdf(fname)

