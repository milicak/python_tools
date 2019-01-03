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

root_folder = '/work/milicak/RUNS/mitgcm/Arctic_4km/'

expid = 'Exp02_0';
dtime = 320;
datadir = root_folder+expid
os.chdir(datadir)
fyear = 1975;

year = 1999;

leapyear = calendar.isleap(year);
yearreg = np.array([31,28,31,30,31,30,31,31,30,31,30,31]);
yearlpy = np.array([31,29,31,30,31,30,31,31,30,31,30,31]);
mnths = ['01','02','03','04','05','06','07','08','09','10','11','12'];

vars = ['ADVr_TH','SALTanom','THETA','UVELMASS','VVELMASS','WVELMASS','WTHMASS']
vars = ['UVELMASS','VVELMASS','WVELMASS','WTHMASS']
# variables = ['ADVr_TH']
# variables = ['SALTanom']
# variables = ['THETA']
# variables = ['UVELMASS']
# variables = ['VVELMASS']
# variables = ['WVELMASS']
# variables = ['WTHMASS']

find0 = 270;
d0 = date(fyear, 1, 1);
d1 = date(year, 1, 1);
delta = d1-d0;

if(leapyear):
    yeardays = yearlpy
else:
    yeardays = yearreg


yeardays = np.array([30,30,30,30,30,30,30,30,30,30,30,30]);
xr_chunks = {'i': 420, 'j': 384}

for variables in vars:
    # itrs = delta.days*find0+np.cumsum(yeardays*find0)
    itrs = (year-fyear)*360*find0+np.cumsum(yeardays*find0)
    print(itrs)
    data = xmitgcm.open_mdsdataset(data_dir=datadir,prefix=variables,iters=itrs,
                                   geometry='curvilinear',read_grid='True',
                                   ignore_unknown_vars='True',
                                   chunks=xr_chunks,
                                   ref_date='1974-12-31 12:0:0',delta_t=dtime)
        
    fname = root_folder + '3DArcticOcean_monthly_'+variables+'_'+np.str(year)+'_1-12.nc'
    # fname = root_folder + '3DArcticOcean_monthly_'+variables[0]+'_'+np.str(year)+'_1-12.nc'
    print(fname)
    print(data)
    data.to_netcdf(fname,engine='netcdf4')
    

os.chdir('/home/ntnu/milicak/python_tools/Analysis/mitgcm/Arctic_4km/Analysis')
