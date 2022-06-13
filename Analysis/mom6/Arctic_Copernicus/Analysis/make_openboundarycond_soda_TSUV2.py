import numpy as np
from os import chdir as cd
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
# # import ESMF
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LinearSegmentedColormap

cd ('/okyanus/users/milicak/dataset/SODA3_12_21/OBC/')

params=[]
params.append({'suffix':'_segment_001','dim0':2,'tr_in':'tracers_north_1_365.nc','tr_out':'obc_ts_north_1_365.nc','uv_in':'uv_north_1_365.nc','uv_out':'obc_uv_north_1_365.nc','ssh_in':'ssh_north_1_365.nc','ssh_out':'obc_ssh_north_1_365.nc'})
params.append({'suffix':'_segment_002','dim0':3,'tr_in':'tracers_west_1_365.nc','tr_out':'obc_ts_west_1_365.nc','uv_in':'uv_west_1_365.nc','uv_out':'obc_uv_west_1_365.nc','ssh_in':'ssh_west_1_365.nc','ssh_out':'obc_ssh_west_1_365.nc'})
params.append({'suffix':'_segment_003','dim0':2,'tr_in':'tracers_south_1_365.nc','tr_out':'obc_ts_south_1_365.nc','uv_in':'uv_south_1_365.nc','uv_out':'obc_uv_south_1_365.nc','ssh_in':'ssh_south_1_365.nc','ssh_out':'obc_ssh_south_1_365.nc'})
params.append({'suffix':'_segment_004','dim0':3,'tr_in':'tracers_east_1_365.nc','tr_out':'obc_ts_east_1_365.nc','uv_in':'uv_east_1_365.nc','uv_out':'obc_uv_east_1_365.nc','ssh_in':'ssh_east_1_365.nc','ssh_out':'obc_ssh_east_1_365.nc'})

for pr in params:
    ds=xr.open_dataset(pr['tr_in'],decode_times=False)
    zl=ds.temp.z_l
    zi=0.5*(np.roll(zl,shift=-1)+zl)
    zi[-1]=6500.
    ds['z_i']=zi
    dz=zi-np.roll(zi,shift=1)
    dz[0]=zi[0]
    ds['dz']=dz
    nt=ds.time.shape[0]
    nx=ds.lon.shape[0]
    dz=np.tile(ds.dz.data[np.newaxis,:,np.newaxis],(nt,1,nx))
    da_dz=xr.DataArray(dz,coords=[('time',ds.time),('z_l',ds.z_l),('locations',ds.locations)])
    da_dz=da_dz.expand_dims('dim_0',pr['dim0'])
    ds.time.attrs['modulo']=' '
    da_temp=xr.DataArray(ds.temp.ffill(dim='locations',limit=None).ffill(dim='z_l').fillna(0.))
    da_temp=da_temp.expand_dims('dim_0',pr['dim0'])
    da_salt=xr.DataArray(ds.salt.ffill(dim='locations',limit=None).ffill(dim='z_l').fillna(0.))
    da_salt=da_salt.expand_dims('dim_0',pr['dim0'])
    ds_=xr.Dataset({'temp'+pr['suffix']:da_temp,'salt'+pr['suffix']:da_salt,'lon':ds.lon,'lat':ds.lat,'dz_temp'+pr['suffix']:da_dz,'dz_salt'+pr['suffix']:da_dz})
    for v in ds_:
        ds_[v].encoding['_FillValue']=1.e20
    ds_.to_netcdf(pr['tr_out'],unlimited_dims=('time'))
    ds=xr.open_dataset(pr['uv_in'],decode_times=False)
    ds.time.attrs['modulo']=' '
    da_u=xr.DataArray(ds.u.ffill(dim='locations',limit=None).ffill(dim='z_l').fillna(0.))
    da_u=da_u.expand_dims('dim_0',pr['dim0'])
    da_u['locations']-=1
    da_v=xr.DataArray(ds.v.ffill(dim='locations',limit=None).ffill(dim='z_l').fillna(0.))
    da_v=da_v.expand_dims('dim_0',pr['dim0'])
    da_v['locations']-=1
    ds_=xr.Dataset({'u'+pr['suffix']:da_u,'v'+pr['suffix']:da_v,'lon':ds.lon,'lat':ds.lat,'dz_u'+pr['suffix']:da_dz,'dz_v'+pr['suffix']:da_dz})
    for v in ds_:
        ds_[v].encoding['_FillValue']=1.e20
    ds_.to_netcdf(pr['uv_out'],unlimited_dims=('time'))
    ds=xr.open_dataset(pr['ssh_in'],decode_times=False)
    ds.time.attrs['modulo']=' '
    da_ssh=xr.DataArray(ds.ssh.ffill(dim='locations',limit=None).fillna(0.))
    da_ssh=da_ssh.expand_dims('dim_0',pr['dim0']-1)
    ds_=xr.Dataset({'ssh'+pr['suffix']:da_ssh,'lon':ds.lon,'lat':ds.lat})
    for v in ds_:
        ds_[v].encoding['_FillValue']=1.e20
    ds_.to_netcdf(pr['ssh_out'],unlimited_dims=('time'))



