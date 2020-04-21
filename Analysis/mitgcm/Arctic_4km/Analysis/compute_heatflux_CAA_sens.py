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
                                                                       
fyear = 1992                                                           
# fyear = 1992                                                         
lyear = 2018                                                           
# lyear = fyear+1                                                      
datadir = root_folder+expid                                            
# os.chdir(datadir)                                                    
prename = '3DArcticOcean_monthly_UTHMASS_' 

fname = datadir + '/' + prename +'*.nc'
list=sorted(glob.glob(fname))

time = pd.date_range('1992-01-01', freq='M', periods=12 * 25)

# df = xr.open_mfdataset(list)
df = xr.open_mfdataset(list,combine='by_coords')     
# old way 
time = pd.date_range('1992-01-01', freq='M', periods=12 * 25)
# new way
# arr = np.array([datetime.datetime(1991, 12, 16) + datetime.timedelta(days=i-1) for i in range(30,365*25+7,30)])
arr = np.array([datetime.datetime(1992, 1, 1) + datetime.timedelta(days=i-1) for i in range(30,365*25+7,30)])
arr = arr[:300]
newtime = xr.Dataset({'time': arr}) 
# df['time'] = newtime['time']

# heat transport at fram strait Cm^3/s
ds = df.UTHMASS*df.dyG*df.drF

ht = ds.data[:,:,1061,788]
HT_fram = ht.sum(axis=(1))
dnm = HT_fram.compute()

# data = xr.DataArray(dnm, dims=('time'), coords={'time': time})
data = xr.DataArray(dnm, dims=('time'), coords={'time': df['time']})
dsnew = data.to_dataset(name='heat_transport')
fname = expid + '_CAA_heat_transport.nc'
dsnew.to_netcdf(fname)

