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
# fyear = 1992
lyear = 2018
# lyear = fyear+1
datadir = root_folder+expid
# os.chdir(datadir)
prename = '2DArcticOcean_'
postname = '_avg.nc'
dropvars = ['SIuice','SIvice','oceTAUY','ETAN','SALT','SIqneti','oceSPDep',
           'SIhsalt','momVort3','oceSPflx','SIheff','THETA','SIqneto','SIhsnow',
           'oceTAUX','dxG','dyG','rAz','dxC','dyC','rAw','rAs']

ens_list = []
SI_extent = np.zeros((lyear-fyear)*12)
k=0
for year in range(fyear,lyear):
    for month in range(1,13):
        sdate = "%2.2d" % (month)
        fname = datadir+'/'+prename+np.str(year)+'_'+sdate+'_avg.nc'
        print(fname)
        ds1 = xr.open_dataset(fname, chunks={'i':500, 'j':500})['SIarea']
                              # drop_variables=dropvars)
        ds1.data[ds1.data>0.15] = 1
        ds1.data[ds1.data<=0.15] = 0
        si = ds1*ds1.rA
        SI_extent[k]=si.sum(dim=['i','j']).compute()
        k+=1
        # ens_list.append(ds1)
        # os.system(cmmnd)



SI_extent = np.reshape(SI_extent,[lyear-fyear,12])
df = pd.DataFrame(SI_extent)
df.to_csv("SI_extent_ctrl_1992_2017")

# ds = xr.concat(ens_list, dim='ensemble')
# os.chdir('/cluster/home/milicak/python_tools/Analysis/mitgcm/Arctic4km/Analysis')
# si = ds.SIarea
# si.data[si.data>0.15] = 1
# si.data[si.data<=0.15] = 0
# dnm = si*ds.rA
# dd = dnm.sum(dim=['i','j'])
# dd.mean().compute()/1e12
