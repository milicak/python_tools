import numpy as np
from os import chdir as cd
import numpy.ma as ma
import glob
import xarray as xr
import sys

year = 2003
root_folder = '/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OBC2/'   

params=[]
params.append({'suffix':'_segment_002','dim0':2,
    'tr_in':'tracers_north','tr_out':root_folder+'obc_ts_north_'+str(year)+'.nc',
    'uv_in':'uv_north','uv_out':root_folder+'obc_uv_north_'+str(year)+'.nc',
    'ssh_in':'ssh_north','ssh_out':root_folder+'obc_ssh_north_'+str(year)+'.nc'})
params.append({'suffix':'_segment_001','dim0':2,
    'tr_in':'tracers_south','tr_out':root_folder+'obc_ts_south_'+str(year)+'.nc',
    'uv_in':'uv_south','uv_out':root_folder+'obc_uv_south_'+str(year)+'.nc',
    'ssh_in':'ssh_south','ssh_out':root_folder+'obc_ssh_south_'+str(year)+'.nc'})

for pr in params:
    print(pr)
    ls0 = sorted(glob.glob(root_folder+str(year-1)+'/'+pr['tr_in'] +'*.nc'))
    ls1 = sorted(glob.glob(root_folder+str(year)+'/'+pr['tr_in'] +'*.nc'))
    ls2 = sorted(glob.glob(root_folder+str(year+1)+'/'+pr['tr_in'] +'*.nc'))
    ls1 = ls0[-1].split() + ls1 + ls2[0].split()
    # ds=xr.open_dataset(pr['tr_in'],decode_times=False)
    ds=xr.open_mfdataset(ls1,decode_times=False)
    ds['locations']=ds.locations.astype('int32')
    zl=ds.temp.depth
    zi=0.5*(np.roll(zl,shift=-1)+zl)
    zi[-1]=6500.
    ds['z_i']=zi
    dz=zi-np.roll(zi,shift=1)
    dz[0]=zi[0]
    ds['dz']=dz
    nt=ds.time.shape[0]
    nx=ds.lon.shape[0]
    dz=np.tile(ds.dz.data[np.newaxis,:,np.newaxis],(nt,1,nx))
    # da_dz=xr.DataArray(dz,coords=[('time',ds.time.data),('z_l',ds.depth.data),('locations',ds.locations.data)])
    da_dz=xr.DataArray(dz,coords=[('time',ds.time.data),('z_l',ds.depth.data),('locations',np.copy(ds.locations.astype('int32')))])
    da_dz=da_dz.expand_dims('dim_0',pr['dim0'])
    # ds.time.attrs['modulo']=' '
    da_temp=xr.DataArray(ds.temp.ffill(dim='locations',limit=None).ffill(dim='depth').fillna(0.))
    da_temp=da_temp.expand_dims('dim_0',pr['dim0'])
    da_salt=xr.DataArray(ds.salt.ffill(dim='locations',limit=None).ffill(dim='depth').fillna(0.))
    da_salt=da_salt.expand_dims('dim_0',pr['dim0'])
    ds_=xr.Dataset({'temp'+pr['suffix']:da_temp,'salt'+pr['suffix']:da_salt,'lon':ds.lon,'lat':ds.lat,'dz_temp'+pr['suffix']:da_dz,'dz_salt'+pr['suffix']:da_dz})
    for v in ds_:
        ds_[v].encoding['_FillValue']=1.e20
    ds_.to_netcdf(pr['tr_out'],unlimited_dims=('time'))
    ls0 = sorted(glob.glob(root_folder+str(year-1)+'/'+pr['uv_in'] +'*.nc'))
    ls1 = sorted(glob.glob(root_folder+str(year)+'/'+pr['uv_in'] +'*.nc'))
    ls2 = sorted(glob.glob(root_folder+str(year+1)+'/'+pr['uv_in'] +'*.nc'))
    ls1 = ls0[-1].split() + ls1 + ls2[0].split()
    # ds=xr.open_dataset(pr['uv_in'],decode_times=False)
    ds=xr.open_mfdataset(ls1,decode_times=False)
    ds['locations']=ds.locations.astype('int32')
    # ds.time.attrs['modulo']=' '
    da_u=xr.DataArray(ds.u.ffill(dim='locations',limit=None).ffill(dim='depth').fillna(0.))
    da_u=da_u.expand_dims('dim_0',pr['dim0'])
    da_u['locations']=da_u['locations']-1
    da_v=xr.DataArray(ds.v.ffill(dim='locations',limit=None).ffill(dim='depth').fillna(0.))
    da_v=da_v.expand_dims('dim_0',pr['dim0'])
    da_v['locations']=da_v['locations']-1
    ds_=xr.Dataset({'u'+pr['suffix']:da_u,'v'+pr['suffix']:da_v,'lon':ds.lon,'lat':ds.lat,'dz_u'+pr['suffix']:da_dz,'dz_v'+pr['suffix']:da_dz})
    for v in ds_:
        ds_[v].encoding['_FillValue']=1.e20
    ds_.to_netcdf(pr['uv_out'],unlimited_dims=('time'))
    ls0 = sorted(glob.glob(root_folder+str(year-1)+'/'+pr['ssh_in'] +'*.nc'))
    ls1 = sorted(glob.glob(root_folder+str(year)+'/'+pr['ssh_in'] +'*.nc'))
    ls2 = sorted(glob.glob(root_folder+str(year+1)+'/'+pr['ssh_in'] +'*.nc'))
    ls1 = ls0[-1].split() + ls1 + ls2[0].split()
    # ds=xr.open_dataset(pr['ssh_in'],decode_times=False)
    ds=xr.open_mfdataset(ls1,decode_times=False)
    # ds.time.attrs['modulo']=' '
    da_ssh=xr.DataArray(ds.ssh.ffill(dim='locations',limit=None).fillna(0.))
    da_ssh=da_ssh.expand_dims('dim_0',pr['dim0']-1)
    ds_=xr.Dataset({'ssh'+pr['suffix']:da_ssh,'lon':ds.lon,'lat':ds.lat})
    for v in ds_:
        ds_[v].encoding['_FillValue']=1.e20
    ds_.to_netcdf(pr['ssh_out'],unlimited_dims=('time'))



