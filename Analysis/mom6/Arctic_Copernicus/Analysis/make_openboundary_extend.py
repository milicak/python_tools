import numpy as np
import xarray as xr
import sys
import pandas as pd

time = pd.date_range("2017-12-24-13", freq="5D", periods=3)
# time = pd.date_range("1980-01-03-13", freq="5D", periods=2777)
root_folder = '/okyanus/users/milicak/dataset/SODA3_12_21/Arctic_Copernicus_OBC/'

era5_dict = {
            'obc_uv_south.nc':'obc_uv_south_ext.nc',
            'obc_ts_south.nc':'obc_ts_south_ext.nc',
            'obc_ssh_south.nc':'obc_ssh_south_ext.nc',
            'obc_uv_north.nc':'obc_uv_north_ext.nc',
            'obc_ts_north.nc':'obc_ts_north_ext.nc',
            'obc_ssh_north.nc':'obc_ssh_north_ext.nc',
            'obc_uv_east.nc':'obc_uv_east_ext.nc',
            'obc_ts_east.nc':'obc_ts_east_ext.nc',
            'obc_ssh_east.nc':'obc_ssh_east_ext.nc',
            'obc_uv_west.nc':'obc_uv_west_ext.nc',
            'obc_ts_west.nc':'obc_ts_west_ext.nc',
            'obc_ssh_west.nc':'obc_ssh_west_ext.nc'
            }

for f, f1 in era5_dict.items():
    print(f)
    fname = root_folder + f
    df = xr.open_dataset(fname)
    ds = df.isel(time=slice(2771,2775))
    ds['time'] = time
    ds.time.attrs['standard_name']='time'
    ds.time.attrs['long_name']='time'
    ds.time.attrs['axis']='T'
    ds.time.attrs['modulo']=' '
    ds.time.attrs['calendar']='gregorian'
    ds2 = xr.concat((df,ds),dim='time')
    # ds2.time.attrs['modulo']=' '
    ds2.time.encoding['calendar']='gregorian'
    ds2.time.encoding.update({'dtype':'float64'})
    outname = root_folder + f1
    ds2.to_netcdf(outname)
