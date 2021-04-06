import xarray as xr
import xesmf as xe
import matplotlib.pyplot as plt
import bottleneck
import numpy as np
import subprocess as sp
import os
import glob
import cartopy.crs as ccrs


def open_grid(path,decode_times=False):
    """Return a grid object containing staggered grid locations"""
    grid={}
    grid['ds']=xr.open_dataset(path,decode_times=False)
    grid['ds']=grid['ds'].drop_dims(['ny','nx'])
    grid['ds']=grid['ds'].drop_vars(['tile'])
    grid['nyp']=grid['ds'].nyp.data[-1]+1
    grid['nxp']=grid['ds'].nxp.data[-1]+1
    nxp=grid['nxp'];nyp=grid['nyp']
    grid['h'] = grid['ds'].isel(nxp=slice(1,nxp+1,2),nyp=slice(1,nyp+1,2))
    #The q grid is not symmetric, but Cu and Cv are
    grid['q'] = grid['ds'].isel(nxp=slice(2,nxp+1,2),nyp=slice(2,nyp+1,2))
    grid['Cu'] = grid['ds'].isel(nxp=slice(0,nxp+1,2),nyp=slice(1,nyp+1,2))
    grid['Cv'] = grid['ds'].isel(nxp=slice(1,nxp+1,2),nyp=slice(0,nyp+1,2))
    return grid


def open_dataset(path,fields,grid):
    ds=xr.open_dataset(path,decode_times=False)

    tracer_list=[];uv_list=[]
    for f in fields:
        for fnam,val in zip(f.keys(),f.values()):
            if val=='h':tracer_list.append(fnam)
            if val=='Cu':uv_list.append(fnam)
            if val=='Cv':uv_list.append(fnam)

    ds_tr = xr.merge([ds, grid['h']])
    ds_u= xr.merge([ds,grid['Cu']])
    ds_v= xr.merge([ds,grid['Cv']])
    return {'ds_tr':ds_tr,'ds_u':ds_u,'ds_v':ds_v,'tracers':tracer_list,'uv':uv_list}


def apply_rotation(u,v,angle_dx,time_slice=slice(0)):
    """Rotate from model space to easterly coordinates"""
    #deg_rad=np.pi/180.
    deg_rad=1.0

    if time_slice is not None:
        t=u.time.isel(time=time_slice)
    else:
        t=u.time

    if time_slice is not None:
        ue=np.cos(angle_dx.data*deg_rad)*u.isel(time=time_slice).data-np.sin(angle_dx.data*deg_rad)*v.isel(time=time_slice).data
        vn=np.sin(angle_dx.data*deg_rad)*u.isel(time=time_slice).data+np.cos(angle_dx.data*deg_rad)*v.isel(time=time_slice).data
    else:
        ue=np.cos(angle_dx.data*deg_rad)*u.data-np.sin(angle_dx.data*deg_rad)*v.data
        vn=np.sin(angle_dx.data*deg_rad)*u.data+np.cos(angle_dx.data*deg_rad)*v.data

    us=ue.shape
    nyp=us[2];nxp=us[3]
    ue=xr.DataArray(ue,coords={'i':np.arange(1,nxp+1),'j':np.arange(1,nyp+1),'time':t,'z_l':u.z_l},dims=('time','z_l','j','i'))
    vn=xr.DataArray(vn,coords={'i':np.arange(1,nxp+1),'j':np.arange(1,nyp+1),'time':t,'z_l':u.z_l},dims=('time','z_l','j','i'))


    return ue,vn


def apply_rotation_transpose(ue,vn,angle_dx,time_slice=slice(0)):
    """Rotate from easterly coordinates to model space"""
    #deg_rad=np.pi/180.
    deg_rad=1.0

    if time_slice is not None:
        t=ue.time.isel(time=time_slice)
    else:
        t=ue.time

    if time_slice is not None:
        u=np.cos(angle_dx.data*deg_rad)*ue.isel(time=time_slice).data+np.sin(angle_dx.data*deg_rad)*vn.isel(time=time_slice).data
        v=-np.sin(angle_dx.data*deg_rad)*ue.isel(time=time_slice).data+np.cos(angle_dx.data*deg_rad)*vn.isel(time=time_slice).data
    else:
        u=np.cos(angle_dx.data*deg_rad)*ue.data+np.sin(angle_dx.data*deg_rad)*vn.data
        v=-np.sin(angle_dx.data*deg_rad)*ue.data+np.cos(angle_dx.data*deg_rad)*vn.data
    nyp=u.shape[2]

    u=xr.DataArray(u,coords={'i':np.arange(1,nyp+1),'time':t,'z_l':ue.z_l},dims=('time','z_l','i'))
    v=xr.DataArray(v,coords={'i':np.arange(1,nyp+1),'time':t,'z_l':ue.z_l},dims=('time','z_l','i'))

    return u,v





path_parent_grid='/okyanus/users/milicak/dataset/MOM6/OM4_025/mosaic.v20170622.unpacked/ocean_hgrid.nc'
parent_grid=open_grid(path_parent_grid)
path_regional_grid='/okyanus/users/milicak/dataset/MOM6/Arctic/ocean_hgrid.nc'
regional_grid=open_grid(path_regional_grid)
regional_grid['ds']=regional_grid['ds'].load()


# plot the angles
regional_grid['ds'].angle_dx.plot()

# load the topo file
dsr_topo=xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic/topog.nc')
dsr_topo = xr.merge([dsr_topo, regional_grid['h']])
#dsr_topo.depth.plot(vmax=-50.,vmin=250.,cmap='gist_gray')
plt.pcolormesh(dsr_topo.depth.data,vmin=-50,vmax=250,cmap=plt.cm.gist_gray)
txt=plt.title('Regional Domain')


