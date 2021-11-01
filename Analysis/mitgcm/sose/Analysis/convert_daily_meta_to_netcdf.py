import xmitgcm
import os
import xarray as xr


# compression options
comp = dict(zlib=True, complevel=5)

root_folder = '/work/opa/mi19918/Projects/mitgcm/sose/'
expid = 'Exp03_0'

folderpath = root_folder + expid + '/'

# itrs = np.arange(288,750000,288)
itrs = np.arange(883354,1768435+345,345)
list = sorted(glob.glob('/work/opa/mi19918/Projects/mitgcm/sose/Exp03_0/*SST*.meta')) 
itrs = np.zeros(len(list))
for ind in range(0,len(list)):
    itrs[ind] = int(list[ind][-15:-5])

# varname = 'SST'
# varname = 'SSS'
# varname = 'ETAN'
# varname = 'SIarea'
varname = 'vort2dsurf'

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
    print(day,mnth,year)
    if dates.year[0] >= 2000:
        if not os.path.isfile(outname):
            encoding = {var: comp for var in df.data_vars}
            if bool(encoding):
                print(outname)
                df.to_netcdf(outname, encoding=encoding)

