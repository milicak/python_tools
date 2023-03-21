import numpy as np

df = xr.open_dataset('/data/opa/bs-mod/upstream_bs-nrt/baseline/inicon/T_utss-sdc2019m01.nc')
ds = df.isel(Y=slice(17,261))
ds.to_netcdf('/work/opa/mi19918/Projects/shared/T_utss-sdc2019m01.nc')


df = xr.open_dataset('/data/opa/bs-mod/upstream_bs-nrt/baseline/inicon/S_utss-sdc2019m01.nc')
ds = df.isel(Y=slice(17,261))
ds.to_netcdf('/work/opa/mi19918/Projects/shared/S_utss-sdc2019m01.nc')

