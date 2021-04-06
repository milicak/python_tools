import xmitgcm
import numpy as np
import os
import xarray as xr


# compression options
comp = dict(zlib=True, complevel=5)

root_folder = '/okyanus/users/milicak/models/MITgcm/Projects/sose/'
expid = 'Exp03_0'

folderpath = root_folder + expid + '/'

# itrs = np.arange(288,750000,288)
# dt1 = 300.0
# itrs = np.arange(288,205920+1,288)
dt1 = 250.0

varname = 'SST'
# varname = 'SSS'
# varname = 'ETAN'
# varname = 'SIarea'
# varname = 'vort2dsurf'

# df1 = xmitgcm.open_mdsdataset(folderpath
#                          ,prefix=varname,read_grid=False)

itrs = np.copy(df1.iter[715::])

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

