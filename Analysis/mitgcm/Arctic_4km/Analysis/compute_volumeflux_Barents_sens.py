import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import os
from datetime import date

root_folder = '/shared/projects/uniklima/globclim/milicak/mitgcm/Arctic_4km/'

expid = 'Exp02_0';
# expid = 'Exp02_1';
# expid = 'Exp02_2';

fyear = 1992
# fyear = 1992
lyear = 2018
# lyear = fyear+1
datadir = root_folder+expid
# os.chdir(datadir)
prename = '3DArcticOcean_monthly_VVELMASS_'

fname = datadir + '/' + prename +'*.nc'
list=sorted(glob.glob(fname))

time = pd.date_range('1992-01-01', freq='M', periods=12 * 25)

df = xr.open_mfdataset(list)
# old way 
time = pd.date_range('1992-01-01', freq='M', periods=12 * 25)
# new way
# arr = np.array([datetime.datetime(1991, 12, 16) + datetime.timedelta(days=i-1) for i in range(30,365*25+7,30)])
arr = np.array([datetime.datetime(1992, 1, 1) + datetime.timedelta(days=i-1) for i in range(30,365*25+7,30)])
arr = arr[:300]
newtime = xr.Dataset({'time': arr}) 
df['time'] = newtime['time']

# volume transport at fram strait m^3/s
# dynew = df.dyC.rename({'j_g': 'j', 'i': 'i_g'})   
ds = df.VVELMASS*df.dxG*df.drF

vt = ds.data[:,:,456,205:531]
VT_barents = vt.sum(axis=(1,2))
dnm = VT_barents.compute()

data = xr.DataArray(dnm, dims=('time'), coords={'time': time})
dsnew = data.to_dataset(name='volume_transport')
fname = expid + '_Barents_volume_transport.nc'
dsnew.to_netcdf(fname)



# VT_fram = []
# for year in range(fyear,lyear):
#     fname = datadir+'/'+prename+np.str(year)+'_1-12.nc'
#     print(fname)
#     df = xr.open_dataset(fname, chunks={'i':500, 'j':500})
#     vt = np.copy(df.UVELMASS[:,:,496:664,580]*df.hFacW[:,496:664,580]
#                   *df.drF*np.copy(df.dyC[496:664,580]))
#     VT_fram = np.append(VT_fram,vt.sum(axis=(1,2)))
#
#
# VT_fram= np.reshape(VT_fram,[lyear-fyear,12])
# df = pd.DataFrame(VT_fram)
# if expid=='Exp02_1':
#     df.to_csv("VT_fram_Atlantic_warm_1992_2007")
# elif expid=='Exp02_2':
#     df.to_csv("VT_fram_Pacific_warm_1992_2007")
# elif expid=='Exp02_0':
#     df.to_csv("VT_fram_ctrl_1992_2007")
#
