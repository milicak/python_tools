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

year = 2000;

leapyear = calendar.isleap(year);
yearreg = np.array([31,28,31,30,31,30,31,31,30,31,30,31]);
yearlpy = np.array([31,29,31,30,31,30,31,31,30,31,30,31]);
mnths = ['01','02','03','04','05','06','07','08','09','10','11','12'];

variables = ['EXFhs','EXFhl','KPPhbl','oceQnet','oceQsw','oceFWflx','oceSflux','PHIBOT']

find0 = 270;

d0 = date(fyear, 1, 1);
d1 = date(year, 1, 1);
delta = d1-d0;
find = (delta.days)*find0+find0

yeardays = np.array([30,30,30,30,30,30,30,30,30,30,30,30]);


for var in variables:
    itrs = (year-fyear)*360*find0+np.cumsum(yeardays*find0)
    print(var)
    print(itrs)
    for ind2 in itrs:
        sdate = "%10.10d" % (ind2)
        lsname ='rm '+var+'.'+sdate+'.data'
        # print(lsname)
        os.system(lsname)
        lsname ='rm '+var+'.'+sdate+'.meta'
        # print(lsname)
        os.system(lsname)
#

    
# os.chdir('/home/ntnu/milicak/python_tools/Analysis/mitgcm/Arctic_4km/Analysis')
