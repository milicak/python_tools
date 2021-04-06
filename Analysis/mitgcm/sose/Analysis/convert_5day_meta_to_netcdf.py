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
varname = 'THETA'

folderpath = root_folder + expid + '/'
df1 = xmitgcm.open_mdsdataset(folderpath
                         ,prefix=varname,read_grid=False)

# dt1 = 300.0
# itrs = np.copy(df1.iter[142])
# itrs = np.copy(df1.iter[:142])
dt1 = 250.0
itrs = np.copy(df1.iter[143::])

# itrs = np.arange(1440,750000,1440)
# itrs = np.arange(524160,750000,1440)
# itrs = np.array([485280])

# varname = 'THETA'
# varname = 'SALT'
# varname = 'UVELMASS'
# varname = 'VVELMASS'
# varname = 'WVELMASS'
#

# varnms = ['THETA', 'SALT', 'UVELMASS', 'VVELMASS', 'WVELMASS']
varnms = ['WVELMASS']
# itr = itrs
# for varname in varnms:
#     df = xmitgcm.open_mdsdataset(folderpath
#                              ,prefix=varname,iters=itr,
#                              ref_date="2004-12-31 12:0:0",delta_t=dt1,
#                              read_grid=False)
#     dates = pd.DatetimeIndex(np.copy(df.time.data))
#     year = str(dates.year[0])
#     mnth = str(dates.month[0])
#     day = str(dates.day[0])
#     outname = (root_folder + 'ncfiles/' + expid + '/' + varname + '_' +
#               year + '_' + mnth.zfill(2) + '_' + day.zfill(2)  + '.nc')
#     if dates.year[0] >= 2000:
#         if not os.path.isfile(outname):
#             encoding = {var: comp for var in df.data_vars}
#             if bool(encoding):
#                 print(outname)
#                 df.to_netcdf(outname, encoding=encoding)

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



