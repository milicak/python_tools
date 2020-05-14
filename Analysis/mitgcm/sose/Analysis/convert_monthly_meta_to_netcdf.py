import xmitgcm
import os
import xarray as xr


# compression options
comp = dict(zlib=True, complevel=5)

root_folder = '/cluster/work/users/milicak/RUNS/mitgcm/sose/'
expid = 'Exp01_0'

folderpath = root_folder + expid + '/'

itrs = np.arange(8784,750000,8784)
# itrs = np.arange(288,5000,288)

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
# varname = 'ADVr_TH'
# varname = 'WTHMASS'
# varname = 'ADVx_TH'
# varname = 'ADVy_TH'
# varname = 'DFxE_TH'
# varname = 'DFyE_TH'
# varname = 'ADVx_SLT'
# varname = 'ADVy_SLT'
# varname = 'DFxE_SLT'
# varname = 'DFyE_SLT'
# varname = 'ADVr_SLT'
# varname = 'EXFevap'
# varname = 'EXFpreci'
# varname = 'EXFempmr'


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

