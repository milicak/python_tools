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

expid = 'Exp01_0';
dtime = 100;
# dtime = 80; # for Exp03_0
datadir = root_folder+expid
os.chdir(datadir)
fyear = 2007;

year = 2007;

leapyear = calendar.isleap(year);
yearreg = np.array([31,28,31,30,31,30,31,31,30,31,30,31]);
yearlpy = np.array([31,29,31,30,31,30,31,31,30,31,30,31]);
mnths = ['01','02','03','04','05','06','07','08','09','10','11','12'];

variables_ocn = ['ETAN','SST','vort2dsurf']
variables = variables_ocn
# variables = ['ETAN']

find0 = 864;
d0 = date(fyear, 1, 1);
d1 = date(year, 1, 1);
delta = d1-d0;

if(leapyear):
    yeardays = yearlpy
else:
    yeardays = yearreg

xr_chunks = {'i': 420, 'j': 384}

for var in variables:
    print(var)
    find = (delta.days)*find0+find0
    for ind in range(0,3):
        lind = find+find0*(yeardays[ind]-1);
        itrs = np.linspace(find,lind,yeardays[ind]);
        print(itrs)
        data = xmitgcm.open_mdsdataset(data_dir=datadir,prefix=var,iters=itrs,
                               geometry='sphericalpolar',read_grid='True',
                               ignore_unknown_vars='True',
                               chunks=xr_chunks,
                               ref_date='2006-12-31 12:0:0',delta_t=dtime)

        fname = root_folder + 'Southern_Ocean_ctrl_daily_ocn_'+var+'_'+np.str(year)+'_'+mnths[ind]+'_01cyc.nc'
        print(fname)
        data.to_netcdf(fname,engine='netcdf4')
        find = lind+find0



os.chdir('/home/ntnu/milicak/python_tools/Analysis/mitgcm/mitgcm_sose/Analysis')

