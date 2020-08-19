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
import dask

root_folder = '/archive/milicak/MITgcm_c65/Projects/Arctic_4km/'

expid = 'Exp02_3';
dtime = 320;
datadir = root_folder+expid
os.chdir(datadir)
fyear = 1992;

year = 1992;
fyear2 = 1992
# fyear2 = 2000
# lyear = fyear2+1;
lyear = 2018

leapyear = calendar.isleap(year);
yearreg = np.array([31,28,31,30,31,30,31,31,30,31,30,31]);
yearlpy = np.array([31,29,31,30,31,30,31,31,30,31,30,31]);
mnths = ['01','02','03','04','05','06','07','08','09','10','11','12'];

vars = ['ADVr_TH','SALT','USLTMASS','VTHMASS',
        'UTHMASS','VVELMASS',
        'UVELMASS','WTHMASS','WVELMASS',
        'THETA','VSLTMASS']

# vars = ['ADVr_TH']
# vars = ['SALT']
# vars = ['USLTMASS']
# vars = ['VTHMASS']
# vars = ['UTHMASS']
# vars = ['VVELMASS']
# vars = ['UVELMASS']
# vars = ['WTHMASS']
# vars = ['WVELMASS']
# vars = ['THETA']
vars = ['VSLTMASS']

find0 = 270;

for year in range(fyear2,lyear):

    d0 = date(fyear, 1, 1);
    d1 = date(year, 1, 1);
    delta = d1-d0;

    if(leapyear):
        yeardays = yearlpy
    else:
        yeardays = yearreg


    if year == 2017:
        yeardays = np.array([30,30,30,30])
    else:
        yeardays = np.array([30,30,30,30,30,30,30,30,30,30,30,30]);

    # xr_chunks = {'i': 420, 'j': 384}
    xr_chunks = {'i_g': 420, 'j': 384, 'k': 10}

    for variables in vars:
        # itrs = delta.days*find0+np.cumsum(yeardays*find0)
        itrs = (year-fyear)*360*find0+np.cumsum(yeardays*find0)
        print(itrs)
        data = xmitgcm.open_mdsdataset(data_dir=datadir,prefix=variables,iters=itrs,
                                       geometry='curvilinear',read_grid='True',
                                       ignore_unknown_vars='True',
                                       # chunks=xr_chunks,
                                       ref_date='1991-12-31 12:0:0',delta_t=dtime)

        fname = root_folder + '3DArcticOcean_monthly_'+variables+'_'+np.str(year)+'_1-12.nc'
        # fname = root_folder + '3DArcticOcean_monthly_'+variables[0]+'_'+np.str(year)+'_1-12.nc'
        print(fname)
        print(data)
        data.load()
        # data.to_netcdf(fname,engine='pynio')
        data.to_netcdf(fname,engine='netcdf4')


os.chdir('/home/milicak/python_tools/Analysis/mitgcm/Arctic_4km/Analysis')

