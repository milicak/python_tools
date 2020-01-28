# ''' This routine converts data to monthly output ```
import calendar
import numpy as np  
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import xmitgcm
import os
from datetime import date

root_folder = '/work/users/mil021/RUNS/mitgcm/Arctic_4km/'

expid = 'Exp02_0';
dtime = 320;
datadir = root_folder+expid
os.chdir(datadir)
fyear = 1975;

year = 1976;

leapyear = calendar.isleap(year);
yearreg = np.array([31,28,31,30,31,30,31,31,30,31,30,31]);
yearlpy = np.array([31,29,31,30,31,30,31,31,30,31,30,31]);
mnths = ['01','02','03','04','05','06','07','08','09','10','11','12'];

variables = ['EXFhs','EXFhl','KPPhbl','oceQnet','oceQsw','oceFWflx','oceSflux','PHIBOT']

find0 = 270;
d0 = date(fyear, 1, 1);
d1 = date(year, 1, 1);
delta = d1-d0;

if(leapyear):
    yeardays = yearlpy
else:
    yeardays = yearreg


itrs = delta.days*find0+np.cumsum(yeardays*find0)
data = xmitgcm.open_mdsdataset(data_dir=datadir,prefix=variables,iters=itrs,
                               geometry='curvilinear',read_grid='True',
                               ignore_unknown_vars='True',
                               ref_date='1974-12-31 12:0:0',delta_t=dtime)
    
fname = root_folder + '2DArcticOcean_monthly_'+np.str(year)+'_1-12.nc'
print fname
data.to_netcdf(fname)

os.chdir('/home/mil021/python_tools/Analysis/mitgcm/Arctic_4km/Analysis')
   

