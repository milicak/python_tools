import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import xmitgcm
import os
from datetime import date


root_folder = '/cluster/work/users/milicak/RUNS/mitgcm/Arctic_4km/'

expid = 'Exp02_0';

fyear = 1992
lyear = 2018
# lyear = fyear+1
datadir = root_folder+expid
os.chdir(datadir)
prename = '2DArcticOcean_'

for year in range(fyear,lyear):
    fname = prename+np.str(year)+'_*_avg.nc'
    list = sorted(glob.glob(fname))
    oname = prename+np.str(year)+'_avg_annual.nc'
    cmmnd = 'cdo ensmean '+fname+' '+oname
    print(cmmnd)
    os.system(cmmnd)



os.chdir('/cluster/home/milicak/python_tools/Analysis/mitgcm/Arctic4km/Analysis')
