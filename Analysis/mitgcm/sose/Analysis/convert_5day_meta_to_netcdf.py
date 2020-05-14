import xmitgcm
import os
import xarray as xr


# compression options
comp = dict(zlib=True, complevel=5)

root_folder = '/cluster/work/users/milicak/RUNS/mitgcm/sose/'
expid = 'Exp01_0'

folderpath = root_folder + expid + '/'

itrs = np.arange(1440,750000,1440)
# itrs = np.arange(288,5000,288)

# varname = 'THETA'
# varname = 'SALT'
# varname = 'UVELMASS'
# varname = 'VVELMASS'
varname = 'WVELMASS'
#
for itr in itrs:
    df = xmitgcm.open_mdsdataset(folderpath
                             ,prefix=varname,iters=itr,
                             ref_date="2004-12-31 12:0:0",delta_t=300,
                             read_grid=False)
    dates = pd.DatetimeIndex(np.copy(df.time.data))
    year = str(dates.year[0])
    mnth = str(dates.month[0])
    day = str(dates.day[0])
    outname = (root_folder + 'ncfiles/' + expid + '/' + varname + '_' +
              year + '_' + mnth.zfill(2) + '_' + day.zfill(2)  + '.nc')
    if not os.path.isfile(outname):
        encoding = {var: comp for var in df.data_vars}
        if bool(encoding):
            print(outname)
            df.to_netcdf(outname, encoding=encoding)

