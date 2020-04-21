import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import os
from datetime import date
import datetime

root_folder = '/shared/projects/uniklima/globclim/milicak/mitgcm/Arctic_4km/'

# expid = 'Exp02_0';
# expid = 'Exp02_1';
expid = 'Exp02_2';

print(expid)

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
dfu = xr.open_mfdataset(list,combine='by_coords',chunks={'j':384,'time':12,'i_g':420})     

fname = datadir + '/' + prenamev +'*.nc'
list=sorted(glob.glob(fname))
dfv = xr.open_mfdataset(list,combine='by_coords',chunks={'j_g':384,'time':12,'i':420})     

aa=np.arange(384,526)      
bb=np.arange(386,386+71+1) 
# just v velocity          
# aa[1::2],bb[:-1]           
# both u and v velocity    
# aa[::2],bb[:-1]            

y = xr.DataArray(bb[:-1], dims='points')
x = xr.DataArray(aa[::2], dims='points')  
x2 = xr.DataArray(aa[1::2], dims='points')  

# volume transport at fram strait m^3/s
# dynew = df.dyC.rename({'j_g': 'j', 'i': 'i_g'})   
ds = dfv.VTHMASS*dfv.dxG*dfv.drF 
ds = ds.to_dataset(name='transport')   
vt1 = ds.transport.isel(i=x,j_g=y)
vt2 = ds.transport.isel(i=x2,j_g=y)

ds1 = dfu.UTHMASS*dfu.dyG*dfu.drF
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

dsnew = xr.merge([VT_barents1, VT_barents2, UT_barents1, UT_barents2])   
fname = expid + '_Barents_heat_transport.nc'
dsnew.to_netcdf(fname)
