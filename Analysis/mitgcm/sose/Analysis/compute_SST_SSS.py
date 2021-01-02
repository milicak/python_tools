import numpy as np
import xarray as xr

root_folder = '/archive2/milicak/mitgcm/sose/'
project_name = 'Exp03_0'
print(project_name)
outname = project_name + '_SST_SSS_mean.nc'

fnames = root_folder+project_name+'/*THETA*.nc'
list = sorted(glob.glob(fnames))
df = xr.open_mfdataset(list)
df = df.sel(time=slice('2007','2011'),k=1)
df = df.mean('time')

fnames = root_folder+project_name+'/*SALT*.nc'
list = sorted(glob.glob(fnames))
ds = xr.open_mfdataset(list)
ds = ds.sel(time=slice('2007','2011'),k=1)
ds = ds.mean('time')

df = df.merge(ds)

comp = dict(zlib=True, complevel=5)
encoding = {var: comp for var in df.data_vars}

fname = root_folder + project_name + '/' + outname
print(fname)
df.to_netcdf(fname, encoding=encoding)


