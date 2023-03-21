import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import matplotlib.pyplot as plt

root_folder = '/archive2/milicak/mitgcm/sose/'
expid1 = 'Exp03_0'
grname = root_folder+expid1+'/grid.nc'
gr = xr.open_dataset(grname)

fnames = root_folder + expid1 + '/SIGMA2*'
list = sorted(glob.glob(fnames))
comp = dict(zlib=True, complevel=5)
# df1 = xr.open_mfdataset(list,concat_dim='time')

# for ind in np.arange(0,len(list)):
for ind in np.arange(146,len(list)):
    print(ind)
    inname = list[ind][-13:]
    fname = root_folder + expid1 + '/' + 'SIGMA2_' + inname
    dft = xr.open_dataset(fname)
    ds = dft.sigma2.where(dft.sigma2>9.7850584)
    dft['sigma2'] = ds
    df = dft.mean('XC')
    outname = root_folder + 'tmp3/' + 'sigma2_' + np.str(ind+1).zfill(3) + '.nc'
    encoding = {var: comp for var in df.data_vars}
    df.to_netcdf(outname)


fnames = root_folder + 'tmp3/sigma2*'
list = sorted(glob.glob(fnames))
tmp = np.zeros((104,1260))
ii = 0;
for ind in np.arange(0,len(list)):
    inname = list[ind]
    df1 = xr.open_dataset(inname)
    tmp = tmp + np.copy(df1.sigma2)
    ii += 1

tmp = tmp/ii
dsn = xr.DataArray(tmp, dims=("Z", "YC"), coords={"Z": df1.Z,
                                                  "YC": df1.YC})
dsn = dsn.to_dataset(name='sigma2')
outname1 = root_folder + expid1 + '/' + expid1 + '_sigma2_mean_time_zonal.nc'
dsn.to_netcdf(outname1)




# time = pd.date_range('2007-01-01', freq='1D', periods=2563)
