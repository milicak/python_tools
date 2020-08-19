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

# compression options
comp = dict(zlib=True, complevel=5)

root_folder = '/archive/milicak/MITgcm_c65/Projects/Arctic_4km/'

expid = 'Exp02_3';
dtime = 320;
datadir = root_folder+expid
os.chdir(datadir)
fyear = 1992;
# fyear2 = 1992
fyear2 = 2000
lyear = 2018 #fyear2+4;

for year in range(fyear2,lyear):

    leapyear = calendar.isleap(year);
    # yearreg = np.array([31,28,31,30,31,30,31,31,30,31,30,31]);
    # yearlpy = np.array([31,29,31,30,31,30,31,31,30,31,30,31]);
    yearreg = np.array([31,28,31,30,31,30,31,31,30,31,30,30]);
    yearlpy = np.array([31,29,31,30,31,30,31,31,30,31,30,30]);
    mnths = ['01','02','03','04','05','06','07','08','09','10','11','12'];

    variables_ocn = ['ETAN','SST','vort2dsurf'];
    variables = variables_ocn

    find0 = 270;
    d0 = date(fyear, 1, 1);
    d1 = date(year, 1, 1);
    delta = d1-d0;
    find = (delta.days)*find0+find0

    if(leapyear):
        yeardays = yearlpy
    else:
        yeardays = yearreg


    if year == 2017:
        yeardays = np.array([30,30,30,30])
    else:
        yeardays = np.array([30,30,30,30,30,30,30,30,30,30,30,30]);


    for ind in range(0,12):
        lind = find+find0*(yeardays[ind]-1);
        # itrs = np.linspace(find,lind,yeardays[ind]);
        itrs = (year-fyear)*360*find0+np.cumsum(yeardays*find0)
        print(itrs)
        data = xmitgcm.open_mdsdataset(data_dir=datadir,prefix=variables,iters=itrs,
                                       geometry='curvilinear',read_grid='False',
                                       ref_date='1991-12-31 12:0:0',delta_t=dtime)

        fname = root_folder + '2DArcticOcean_'+np.str(year)+'_'+mnths[ind]+'.nc'
        print(fname)
        encoding = {var: comp for var in data.data_vars}
        data.load()
        data.to_netcdf(fname, encoding=encoding)
        find = lind+find0




os.chdir('/home/milicak/python_tools/Analysis/mitgcm/Arctic_4km/Analysis')

