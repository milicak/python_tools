import numpy as np
import xarray as xr
import os

root_folder = '/archive/milicak/MITgcm_c65/Projects/Arctic_4km/'

expid  = 'Exp02_3'
tmpid = 'tmp'

prefname1 = '2DArcticOcean_'

comp = dict(zlib=True, complevel=5)

for year in range(1992,2018):
    for month in range(1,13):
        print(year,month)
        fname = root_folder + '/' + prefname1 + str(year) + '_' + str(month).zfill(2) + '.nc'
        fname2 = root_folder + tmpid + '/' + prefname1 + str(year) + '_' + str(month).zfill(2) + '.nc'
        df = xr.open_dataset(fname)
        df2 = df.drop(('rA','rAz','XC','YC','XG','YG','CS','SN','Z','Zp1','Zu','Zl','dxG',
                     'dyG','Depth','dxC','dyC','rAw','rAs','drC','drF','PHrefC','PHrefF','hFacC','hFacW','hFacS'))
        encoding = {var: comp for var in df2.data_vars}
        df2.to_netcdf('dnm.nc', encoding=encoding)
        cmnd = 'mv dnm.nc ' + fname2
        os.system(cmnd)

