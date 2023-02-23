'''This is code is written using
   https://github.com/ESMG/PyCNAL_regridding/blob/master/examples/SODA3.3.1/Creating_Initial_and_Boundary_conditions_from_SODA3.py'''
import os
import numpy as np
import xesmf
import glob
import xarray as xr
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings('ignore')

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

def apply_rotation(u,v,angle_dx,time_slice=slice(0)):
    """Rotate from model space to easterly coordinates"""
    deg_rad=np.pi/180.

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
    ue=xr.DataArray(ue,coords={'i':np.arange(1,nxp+1),'j':np.arange(1,nyp+1),'time':t,'depth':u.depth},dims=('time','depth','j','i'))
    vn=xr.DataArray(vn,coords={'i':np.arange(1,nxp+1),'j':np.arange(1,nyp+1),'time':t,'depth':u.depth},dims=('time','depth','j','i'))


    return ue,vn

def apply_rotation_transpose(ue,vn,angle_dx,time_slice=slice(0)):
    """Rotate from easterly coordinates to model space"""
    deg_rad=np.pi/180.

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

    u=xr.DataArray(u,coords={'i':np.arange(1,nyp+1),'time':t,'depth':ue.depth},dims=('time','depth','i'))
    v=xr.DataArray(v,coords={'i':np.arange(1,nyp+1),'time':t,'depth':ue.depth},dims=('time','depth','i'))

    return u,v

def velocity_at_corners(ds_u,ds_v):
    #nxp and nyp are 1 points more than the lon, lat lengths, respectively
    nxp = ds_u.shape[3]+1 ;  nyp = ds_u.shape[2]+1
    #upper-right q points
    u_q=0.5*(ds_u+ds_u.roll(roll_coords='lat',lat=-1)).isel(lon=slice(1,nxp))
    #upper-right q points
    v_q=0.5*(ds_v+ds_v.roll(roll_coords='lon',lon=-1)).isel(lat=slice(1,nyp))
    ds_uvq = xr.Dataset({'u':u_q,'v':v_q},coords={'time':ds_u.time,'lon':parent_grid['q'].x,'lat':parent_grid['q'].y,'angle_dx':parent_grid['q'].angle_dx})
    return ds_uvq

root_folder = '/Volumes/A1/workdir/milicak/datasets/MOM6/NWA25/OBC2/'
# root_folder = '/Volumes/A1/workdir/milicak/datasets/MOM6/NWA25/OBC/'
mom_dir = '/Volumes/A1/workdir/milicak/datasets/MOM6/NWA25/'
path_regional_grid = mom_dir + './ocean_hgrid.nc'
regional_grid=open_grid(path_regional_grid)


# Define boundaries for regional domain model
ds_regional=regional_grid['ds']
# northern boundary
north = xr.Dataset()
north['lon'] = ds_regional['x'].isel(nyp=-1)
north['lat'] = ds_regional['y'].isel(nyp=-1)
# southern boundary
south = xr.Dataset()
south['lon'] = ds_regional['x'].isel(nyp=0)
south['lat'] = ds_regional['y'].isel(nyp=0)
# eastern boundary                           
east = xr.Dataset()                          
east['lon'] = ds_regional['x'].isel(nxp=-1)  
east['lat'] = ds_regional['y'].isel(nxp=-1)  


# get an initial dataset to create the weights
year = 1997
ls1 = sorted(glob.glob('/Volumes/A1/workdir/milicak/datasets/glorys/*'+str(year)+'*.nc'))
dfm = xr.open_dataset(ls1[0])
dft = dfm.get(['thetao', 'so', 'zos'])
dft = dft.rename({'longitude':'lon','latitude':'lat'})
# get velocities
dfuv = dfm.get(['uo','vo'])
dfuv = dfuv.rename({'longitude':'lon','latitude':'lat'})


# Calculate remapping weights
# Using nearest neighbor - other options could be used here , e.g. bilinear.
method = 'bilinear'
# method = 'nearest_s2d'
regrid_north_uv = xesmf.Regridder(dfuv, north, method,
                       locstream_out=True, periodic=False, filename='regrid_north_uv.nc',reuse_weights=True)

regrid_south_uv = xesmf.Regridder(dfuv, south, method,
                       locstream_out=True, periodic=False, filename='regrid_south_uv.nc',reuse_weights=True)

regrid_east_uv = xesmf.Regridder(dfuv, east, method,
                       locstream_out=True, periodic=False, filename='regrid_east_uv.nc',reuse_weights=True)

# tracer values
regrid_north_tr = xesmf.Regridder(dft, north, method,
                       locstream_out=True, periodic=False, filename='regrid_north_tr.nc',reuse_weights=True)

regrid_south_tr = xesmf.Regridder(dft, south, method,
                       locstream_out=True, periodic=False, filename='regrid_south_tr.nc',reuse_weights=True)

regrid_east_tr = xesmf.Regridder(dft, east, method,
                       locstream_out=True, periodic=False, filename='regrid_east_tr.nc',reuse_weights=True)

for idx, fname in enumerate(ls1):
    print(idx,fname)
    dfm = xr.open_dataset(fname)
    dft = dfm.get(['thetao', 'so', 'zos'])
    dfuv = dfm.get(['uo','vo'])
    # interpolate T/S and write to netcdf file
    temp_north = regrid_north_tr(dft['thetao'])
    salt_north = regrid_north_tr(dft['so'])
    temp_north = temp_north.bfill(dim='locations').ffill(dim='depth')
    salt_north = salt_north.bfill(dim='locations').ffill(dim='depth')
    # this is necesarry for this particular domain where north boundary has
    # land points at the beginning
    temp_north = temp_north.bfill(dim='locations')
    salt_north = salt_north.bfill(dim='locations')
    ds_tr_north = xr.Dataset({'temp':temp_north,'salt':salt_north})
    ds_tr_north.time.encoding['calendar']='gregorian'
    fnam=root_folder + str(year) + '/' + 'tracers_north_' + np.str(idx).zfill(4)+ '_obc.nc'
    ds_tr_north.to_netcdf(fnam,unlimited_dims='time',format='NETCDF4')

    temp_south = regrid_south_tr(dft['thetao'])
    salt_south = regrid_south_tr(dft['so'])
    temp_south = temp_south.ffill(dim='locations').ffill(dim='depth')
    salt_south = salt_south.ffill(dim='locations').ffill(dim='depth')
    # this is NOT necesarry for this domain for the south boundary
    # I decided to keep to be consistent with north
    temp_south = temp_south.bfill(dim='locations')
    salt_south = salt_south.bfill(dim='locations')
    ds_tr_south = xr.Dataset({'temp':temp_south,'salt':salt_south})
    ds_tr_south.time.encoding['calendar']='gregorian'
    fnam=root_folder + str(year) + '/' + 'tracers_south_' + np.str(idx).zfill(4)+ '_obc.nc'
    ds_tr_south.to_netcdf(fnam,unlimited_dims='time',format='NETCDF4')

    temp_east = regrid_east_tr(dft['thetao'])
    salt_east = regrid_east_tr(dft['so'])
    temp_east = temp_east.ffill(dim='locations').ffill(dim='depth')
    salt_east = salt_east.ffill(dim='locations').ffill(dim='depth')
    # this is NOT necesarry for this domain for the east boundary
    # I decided to keep to be consistent with north
    temp_east = temp_east.bfill(dim='locations')
    salt_east = salt_east.bfill(dim='locations')
    ds_tr_east = xr.Dataset({'temp':temp_east,'salt':salt_east})
    ds_tr_east.time.encoding['calendar']='gregorian'
    fnam=root_folder + str(year) + '/' + 'tracers_east_' + np.str(idx).zfill(4)+ '_obc.nc'
    ds_tr_east.to_netcdf(fnam,unlimited_dims='time',format='NETCDF4')

    # interpolate ssh
    ssh_north = regrid_north_tr(dft['zos'])
    ssh_north = ssh_north.bfill(dim='locations').ffill(dim='locations')
    ds_ssh_north = xr.Dataset({'ssh':ssh_north})
    ds_ssh_north.time.encoding['calendar'] = "gregorian"
    fnam=root_folder + str(year) + '/' + 'ssh_north' + np.str(idx).zfill(4)+ '_obc.nc'
    ds_ssh_north.to_netcdf(fnam,unlimited_dims='time',format='NETCDF4')

    ssh_south = regrid_south_tr(dft['zos'])
    ssh_south = ssh_south.bfill(dim='locations').ffill(dim='locations')
    # ssh_south = ssh_south.ffill(dim='locations').bfill(dim='locations')
    ds_ssh_south = xr.Dataset({'ssh':ssh_south})
    ds_ssh_south.time.encoding['calendar'] = "gregorian"
    fnam=root_folder+ str(year) + '/' + 'ssh_south' + np.str(idx).zfill(4)+ '_obc.nc'
    ds_ssh_south.to_netcdf(fnam,unlimited_dims='time',format='NETCDF4')

    ssh_east = regrid_east_tr(dft['zos'])
    ssh_east = ssh_east.ffill(dim='locations').bfill(dim='locations')
    ds_ssh_east = xr.Dataset({'ssh':ssh_east})
    ds_ssh_east.time.encoding['calendar'] = "gregorian"
    fnam=root_folder+ str(year) + '/' + 'ssh_east' + np.str(idx).zfill(4)+ '_obc.nc'
    ds_ssh_east.to_netcdf(fnam,unlimited_dims='time',format='NETCDF4')

    # interpolate true u and v velocities
    u_north_r = regrid_north_uv(dfuv['uo'])
    v_north_r = regrid_north_uv(dfuv['vo'])
    u_south_r = regrid_south_uv(dfuv['uo'])
    v_south_r = regrid_south_uv(dfuv['vo'])
    u_east_r = regrid_east_uv(dfuv['uo'])
    v_east_r = regrid_east_uv(dfuv['vo'])
    # rotate back to model grid
    u_north,v_north=apply_rotation_transpose(u_north_r,v_north_r,ds_regional.angle_dx.isel(nyp=ds_regional.nyp[-1]),time_slice=None)
    u_north = u_north.fillna(0)
    v_north = v_north.fillna(0)
    ds_uv_north = xr.Dataset({'u':u_north,'v':v_north},coords={'lon':north.lon,'lat':north.lat,'depth':dfuv.depth})
    ds_uv_north = ds_uv_north.rename({'i':'locations'})
    ds_uv_north.time.encoding['calendar']='gregorian'
    fnam=root_folder + str(year) + '/' + 'uv_north' + np.str(idx).zfill(4)+ '_obc.nc'
    ds_uv_north.to_netcdf(fnam,unlimited_dims='time',format='NETCDF4')

    u_south,v_south=apply_rotation_transpose(u_south_r,v_south_r,ds_regional.angle_dx.isel(nyp=ds_regional.nyp[0]),time_slice=None)
    u_south = u_south.fillna(0)
    v_south = v_south.fillna(0)
    ds_uv_south = xr.Dataset({'u':u_south,'v':v_south},coords={'lon':south.lon,'lat':south.lat,'depth':dfuv.depth})
    ds_uv_south = ds_uv_south.rename({'i':'locations'})
    ds_uv_south.time.encoding['calendar']='gregorian'
    fnam=root_folder + str(year) + '/' + 'uv_south' + np.str(idx).zfill(4)+ '_obc.nc'
    ds_uv_south.to_netcdf(fnam,unlimited_dims='time',format='NETCDF4')

    u_east,v_east=apply_rotation_transpose(u_east_r,v_east_r,ds_regional.angle_dx.isel(nxp=ds_regional.nxp[-1]),time_slice=None)
    u_east = u_east.fillna(0)
    v_east = v_east.fillna(0)
    ds_uv_east = xr.Dataset({'u':u_east,'v':v_east},coords={'lon':east.lon,'lat':east.lat,'depth':dfuv.depth})
    ds_uv_east = ds_uv_east.rename({'i':'locations'})
    ds_uv_east.time.encoding['calendar']='gregorian'
    fnam=root_folder + str(year) + '/' + 'uv_east' + np.str(idx).zfill(4)+ '_obc.nc'
    ds_uv_east.to_netcdf(fnam,unlimited_dims='time',format='NETCDF4')




