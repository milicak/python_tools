import xmitgcm
import os
import xarray as xr


# compression options
comp = dict(zlib=True, complevel=5)

root_folder = '/cluster/work/users/milicak/RUNS/mitgcm/sose/'
expid = 'Exp01_0'

folderpath = root_folder + expid + '/'

# itrs = np.arange(288,750000,288)
# itrs = np.arange(525312,750000,288)
# itrs = 288
itrs = np.array([736128])

varname = 'SST'
# varname = 'SSS'
# varname = 'ETAN'
# varname = 'SIarea'
# varname = 'vort2dsurf'

for itr in itrs:
    df = xmitgcm.open_mdsdataset(folderpath
                             ,prefix=varname,iters=itr,
                             ref_date="1998-01-02 12:0:0",delta_t=300,
                             read_grid=True)
    dates = pd.DatetimeIndex(np.copy(df.time.data))
    year = str(dates.year[0])
    mnth = str(dates.month[0])
    day = str(dates.day[0])
    outname = (root_folder + 'ncfiles/' + expid + '/grid.nc')
    encoding = {var: comp for var in df.data_vars}
    df.to_netcdf(outname, encoding=encoding)

