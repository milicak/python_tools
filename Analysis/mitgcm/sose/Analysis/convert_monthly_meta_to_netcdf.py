import numpy as np
import xmitgcm
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import os
import xarray as xr


# compression options
comp = dict(zlib=True, complevel=5)

root_folder = '/work/opa/mi19918/Projects/mitgcm/sose/'
expid = 'Exp03_0'

folderpath = root_folder + expid + '/'

# itrs = np.arange(8784,750000,8784)
# itrs = np.arange(518256,750000,8784)
list = sorted(glob.glob('/work/opa/mi19918/Projects/mitgcm/sose/Exp03_0/*oceTAUX*.meta')) 
itrs = np.zeros(len(list))
for ind in range(0,len(list)):
    itrs[ind] = int(list[ind][-15:-5])

# varname = 'oceTAUX'
# varname = 'oceTAUY'
# varname = 'oceFWflx'
# varname = 'oceSflux'
# varname = 'oceQnet'
# varname = 'SIheff'
# varname = 'SIuice'
# varname = 'SIvice'
# varname = 'SIhsnow'
# varname = 'SIqneti'
# varname = 'SIqneto'
# varname = 'EXFhl'
# varname = 'EXFhs'
# varname = 'EXFevap'
# varname = 'EXFpreci'
# varname = 'EXFempmr'

# varname = 'WTHMASS'
# varname = 'ADVx_TH'
# varname = 'ADVy_TH'
# varname = 'ADVr_TH'
# varname = 'DFxE_TH'
# varname = 'DFyE_TH'
# varname = 'ADVx_SLT'
# varname = 'ADVy_SLT'
# varname = 'DFxE_SLT'
# varname = 'DFyE_SLT'
varname = 'ADVr_SLT'


for itr in itrs:
    df = xmitgcm.open_mdsdataset(folderpath
                             ,prefix=varname,iters=itr,
                             ref_date="1998-01-02 12:0:0",delta_t=250,
                             read_grid=False,nx=4320,ny=1260,nz=104)
    dates = pd.DatetimeIndex(np.copy(df.time.data))
    year = str(dates.year[0])
    mnth = str(dates.month[0])
    day = str(dates.day[0])
    outname = (root_folder + 'ncfiles/' + expid + '/' + varname + '_2nd_' +
              year + '_' + mnth.zfill(2) + '_' + day.zfill(2)  + '.nc')
    if dates.year[0] >= 2000:
        if not os.path.isfile(outname):
            encoding = {var: comp for var in df.data_vars}
            if bool(encoding):
                print(outname)
                df.to_netcdf(outname, encoding=encoding)


