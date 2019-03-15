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

root_folder = '/cluster/work/users//milicak/RUNS/mitgcm/Arctic_4km/Exp02_0/'

datadir = root_folder
os.chdir(datadir)
fyear = 1992;
lyear = 2018;

prename = '3DArcticOcean_monthly_THETA_'

for year in range(fyear,lyear):
    fname = prename+np.str(year)+'_1-12.nc'
    cmmnd = 'ncap2 -s "where(THETA < -2.0) THETA=-2.0" '+fname+ ' tmp/dnm.nc'
    print(cmmnd)
    os.system(cmmnd)
    cmmnd = 'mv tmp/dnm.nc tmp/'+fname
    print(cmmnd)
    os.system(cmmnd)


os.chdir('/cluster/home/milicak/python_tools/Analysis/mitgcm/Arctic4km/Analysis')
