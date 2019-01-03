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

datadir = root_folder
os.chdir(datadir)
fyear = 1975;
lyear = 2000;
skipyear = 17;

prename = '2DArcticOcean_monthly_'

for year in range(fyear,lyear):
    fname = prename+np.str(year)+'_1-12.nc'
    oname = prename+np.str(year+skipyear)+'_1-12.nc'
    cmmnd = 'mv '+fname+ ' dnm/'+oname
    print(cmmnd)
    os.system(cmmnd)


os.chdir('/home/ntnu/milicak/python_tools/Analysis/mitgcm/Arctic_4km/Analysis')
