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

root_folder = '/work/milicak/RUNS/mitgcm/mitgcm_sose/'

expid = 'Exp03_0';
dtime = 80;
# dtime = 100;
datadir = root_folder+expid
os.chdir(datadir)
fyear = 2007;

year = 2014;

leapyear = calendar.isleap(year);
yearreg = np.array([31,28,31,30,31,30,31,31,30,31,30,31]);
yearlpy = np.array([31,29,31,30,31,30,31,31,30,31,30,31]);
mnths = ['01','02','03','04','05','06','07','08','09','10','11','12'];

variables_ocn = ['oceTAUX','oceTAUY','oceFWflx','oceSflux','oceQnet','KPPhbl','EXFhl',
                'EXFhs','SALT','THETA','UVELMASS','VVELMASS','WVELMASS','ADVr_TH','WTHMASS']
variables_ice =['SIarea','SIheff','SIuice','SIvice','SIhsnow','SIqneti','SIqneto']
variables = variables_ocn + variables_ice
# variables = variables_ice
# variables = ['SALT']
# variables = ['oceTAUX']

find0 = 86400/dtime;
d0 = date(fyear, 1, 1);
d1 = date(year, 1, 1);
delta = d1-d0;

if(leapyear):
    yeardays = yearlpy
else:
    yeardays = yearreg


# yeardays = np.array([30,30,30,30,30,30,30,30,30,30,30,30]);
xr_chunks = {'i': 432, 'j': 320}

for var in variables:
    print(var)
    # itrs = (year-fyear)*360*find0+np.cumsum(yeardays*find0)
    itrs = delta.days*find0+np.cumsum(yeardays*find0)
    print(itrs)
    data = xmitgcm.open_mdsdataset(data_dir=datadir,prefix=var,iters=itrs,
                               geometry='sphericalpolar',read_grid='True',
                               ignore_unknown_vars='True',
                               # chunks=xr_chunks,
                               ref_date='2006-12-31 12:0:0',delta_t=dtime)

    fname = root_folder + 'Southern_Ocean_ctrl_monthly_ocn_'+var+'_'+np.str(year)+'_1-12_01cyc.nc'
    data.load()
    print(fname)
    data.to_netcdf(fname,engine='netcdf4')



os.chdir('/home/ntnu/milicak/python_tools/Analysis/mitgcm/mitgcm_sose/Analysis')

