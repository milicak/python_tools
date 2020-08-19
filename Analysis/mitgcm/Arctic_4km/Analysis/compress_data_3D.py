import numpy as np
import os

root_folder = '/archive/milicak/MITgcm_c65/Projects/Arctic_4km/'

expid  = 'Exp02_3'
tmpid = 'tmp'

prefname1 = '3DArcticOcean_monthly_'

vars = ['ADVr_TH','SALT','USLTMASS','VTHMASS',
        'UTHMASS','VVELMASS',
        'UVELMASS','WTHMASS','WVELMASS',
        'THETA','VSLTMASS']

comp = dict(zlib=True, complevel=5)


# for year in range(1992,2018):
for year in range(2009,2018):
    for variables in vars:
        print(year,variables)
        fname = root_folder + '/' + prefname1 + variables + '_' + str(year) + '_1-12.nc'
        fname2 = root_folder + tmpid + '/' + prefname1 + variables + '_' + str(year) + '_1-12.nc'
        outname = root_folder + tmpid + '/dnm2.nc'
        df = xr.open_dataset(fname)
        df2 = df.drop(('rA','rAz','XC','YC','XG','YG','CS','SN','Z','Zp1','Zu','Zl','dxG',
                     'dyG','Depth','dxC','dyC','rAw','rAs','drC','drF','PHrefC','PHrefF','hFacC','hFacW','hFacS'))
        encoding = {var: comp for var in df2.data_vars}
        df2.to_netcdf(outname, encoding=encoding)
        cmnd = 'mv '+ outname + ' ' + fname2
        os.system(cmnd)

