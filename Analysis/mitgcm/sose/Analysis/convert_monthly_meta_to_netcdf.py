import xmitgcm
import numpy as np
import numpy.ma as ma
import glob
import sys
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LinearSegmentedColormap
import os
import xarray as xr


# compression options
comp = dict(zlib=True, complevel=5)

root_folder = '/okyanus/users/milicak/models/MITgcm/Projects/sose/'
expid = 'Exp03_0'

varname = 'oceTAUX'
folderpath = root_folder + expid + '/'
df1 = xmitgcm.open_mdsdataset(folderpath
                         ,prefix=varname,read_grid=False)

# dt1 = 300.0
# itrs = np.copy(df1.iter[:23])
dt1 = 250.0
itrs = np.copy(df1.iter[23::])
#
varname = 'oceTAUX'
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
# varname = 'ADVr_SLT'

varnms = ['oceTAUX', 'oceTAUY', 'oceFWflx', 'oceSflux', 'oceQnet', 'SIheff',
          'SIuice', 'SIvice', 'SIhsnow', 'SIqneti', 'SIqneto', 'EXFhl',
          'EXFhs', 'EXFevap', 'EXFpreci', 'EXFempmr', 'WTHMASS', 'ADVx_TH',
          'ADVy_TH', 'ADVr_TH', 'DFxE_TH', 'DFyE_TH', 'ADVx_SLT', 'ADVy_SLT',
          'DFxE_SLT', 'DFyE_SLT', 'DFyE_SLT']

for varname in varnms:
    for itr in itrs:
        df = xmitgcm.open_mdsdataset(folderpath
                                 ,prefix=varname,iters=itr,
                                 ref_date="2004-12-31 12:0:0",delta_t=dt1,
                                 read_grid=False)
        dates = pd.DatetimeIndex(np.copy(df.time.data))
        year = str(dates.year[0])
        mnth = str(dates.month[0])
        day = str(dates.day[0])
        outname = (root_folder + 'ncfiles/' + expid + '/' + varname + '_' +
                  year + '_' + mnth.zfill(2) + '_' + day.zfill(2)  + '.nc')
        if dates.year[0] >= 2000:
            if not os.path.isfile(outname):
                encoding = {var: comp for var in df.data_vars}
                if bool(encoding):
                    print(outname)
                    df.to_netcdf(outname, encoding=encoding)


