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

variables_ocn = ['ETAN','oceSPDep','oceSPflx','oceTAUX','oceTAUY','SSS','SST','vort2dsurf'];
variables_ice = ['SIarea','SIheff','SIhsalt','SIhsnow','SIqneti','SIqneto','SIuice','SIvice'];
variables = variables_ocn + variables_ice

find0 = 270;
d0 = date(fyear, 1, 1);
d1 = date(year, 1, 1);
delta = d1-d0;
find = (delta.days)*find0+find0

if(leapyear):
    yeardays = yearlpy
else:
    yeardays = yearreg


for ind in xrange(0,12):
    lind = find+find0*(yeardays[ind]-1);
    itrs = np.linspace(find,lind,yeardays[ind]);
    data = xmitgcm.open_mdsdataset(data_dir=datadir,prefix=variables,iters=itrs,
                                   geometry='curvilinear',read_grid='True',
                                   ref_date='1974-12-31 12:0:0',delta_t=dtime)
    
    fname = root_folder + '2DArcticOcean_'+np.str(year)+'_'+mnths[ind]+'.nc'
    print fname
    data.to_netcdf(fname)
    find = lind+find0
    

os.chdir('/home/mil021/python_tools/Analysis/mitgcm/Arctic_4km/Analysis')
