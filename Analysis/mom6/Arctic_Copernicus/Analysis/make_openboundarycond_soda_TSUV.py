'''This is code is written using
   https://github.com/ESMG/PyCNAL_regridding/blob/master/examples/SODA3.3.1/Creating_Initial_and_Boundary_conditions_from_SODA3.py'''
import os
import numpy as np
import xesmf as xe
import xesmf
import glob
import xarray as xr
import scipy.io
from scipy.io import savemat
from scipy.io import loadmat
# from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.colors as colors
from scipy.signal import medfilt2d
import netCDF4
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.path import Path
#for interpolation
from scipy.spatial import cKDTree
from HCtFlood.kara import flood_kara
# from PyCNAL_regridding import *


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

    u=xr.DataArray(u,coords={'i':np.arange(1,nyp+1),'time':t,'st_ocean':ue.st_ocean},dims=('time','st_ocean','i'))
    v=xr.DataArray(v,coords={'i':np.arange(1,nyp+1),'time':t,'st_ocean':ue.st_ocean},dims=('time','st_ocean','i'))

    return u,v


root_folder = '/okyanus/users/milicak/dataset/SODA3.12.21/'
fname1 = 'soda3.12.2_5dy_ocean_reg_1980_01_03.nc'
df = xr.open_dataset(root_folder + fname1)
mom_dir = '/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/'
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
# western boundary
west = xr.Dataset()
west['lon'] = ds_regional['x'].isel(nxp=0)
west['lat'] = ds_regional['y'].isel(nxp=0)
# eastern boundary
east = xr.Dataset()
east['lon'] = ds_regional['x'].isel(nxp=-1)
east['lat'] = ds_regional['y'].isel(nxp=-1)


ls1 = sorted(glob.glob('/okyanus/users/milicak/dataset/SODA3.12.21/*.nc'))
dfm = xr.open_dataset(ls1[0])
dft = dfm.get(['temp', 'salt', 'ssh'])
dfuv = dfm.get(['u','v'])
dft = dft.rename({'xt_ocean':'lon','yt_ocean':'lat'})
dfuv = dfuv.rename({'xu_ocean':'lon','yu_ocean':'lat'})


# Calculate remapping weights
# Using nearest neighbor - other options could be used here , e.g. bilinear.
regrid_north_uv = xesmf.Regridder(dfuv, north, 'nearest_s2d',
                       locstream_out=True, periodic=False, filename='regrid_north_uv.nc',reuse_weights=True)

regrid_south_uv = xesmf.Regridder(dfuv, south, 'nearest_s2d',
                       locstream_out=True, periodic=False, filename='regrid_south_uv.nc',reuse_weights=True)

regrid_east_uv = xesmf.Regridder(dfuv, east, 'nearest_s2d',
                       locstream_out=True, periodic=False, filename='regrid_east_uv.nc',reuse_weights=True)

regrid_west_uv = xesmf.Regridder(dfuv, west, 'nearest_s2d',
                       locstream_out=True, periodic=False, filename='regrid_west_uv.nc',reuse_weights=True)

# tracer values
regrid_north_tr = xesmf.Regridder(dft, north, 'nearest_s2d',
                       locstream_out=True, periodic=False, filename='regrid_north_tr.nc',reuse_weights=True)

regrid_south_tr = xesmf.Regridder(dft, south, 'nearest_s2d',
                       locstream_out=True, periodic=False, filename='regrid_south_tr.nc',reuse_weights=True)

regrid_east_tr = xesmf.Regridder(dft, east, 'nearest_s2d',
                       locstream_out=True, periodic=False, filename='regrid_east_tr.nc',reuse_weights=True)

regrid_west_tr = xesmf.Regridder(dft, west, 'nearest_s2d',
                       locstream_out=True, periodic=False, filename='regrid_west_tr.nc',reuse_weights=True)


for idx, fname in enumerate(ls1):
    print(idx,fname)
    dfm = xr.open_dataset(fname)
    dft = dfm.get(['temp', 'salt', 'ssh'])
    dfuv = dfm.get(['u','v'])
    # interpolate T/S and write to netcdf file
    temp_north = regrid_north_tr(dft['temp'])
    salt_north = regrid_north_tr(dft['salt'])
    temp_north = temp_north.ffill(dim='locations').ffill(dim='st_ocean')
    salt_north = salt_north.ffill(dim='locations').ffill(dim='st_ocean')
    ds_tr_north = xr.Dataset({'temp':temp_north,'salt':salt_north})
    ds_tr_north.time.encoding['calendar']='gregorian'
    fnam=root_folder + 'tracers_north_' + np.str(idx).zfill(4)+ '_obc.nc'
    ds_tr_north.to_netcdf(fnam,unlimited_dims='time',format='NETCDF3_CLASSIC')

    temp_south = regrid_south_tr(dft['temp'])
    salt_south = regrid_south_tr(dft['salt'])
    temp_south = temp_south.ffill(dim='locations').ffill(dim='st_ocean')
    salt_south = salt_south.ffill(dim='locations').ffill(dim='st_ocean')
    ds_tr_south = xr.Dataset({'temp':temp_south,'salt':salt_south})
    ds_tr_south.time.encoding['calendar']='gregorian'
    fnam=root_folder + 'tracers_south_' + np.str(idx).zfill(4)+ '_obc.nc'
    ds_tr_south.to_netcdf(fnam,unlimited_dims='time',format='NETCDF3_CLASSIC')

    temp_east = regrid_east_tr(dft['temp'])
    salt_east = regrid_east_tr(dft['salt'])
    temp_east = temp_east.ffill(dim='locations').ffill(dim='st_ocean')
    salt_east = salt_east.ffill(dim='locations').ffill(dim='st_ocean')
    ds_tr_east = xr.Dataset({'temp':temp_east,'salt':salt_east})
    ds_tr_east.time.encoding['calendar']='gregorian'
    fnam=root_folder + 'tracers_east' + np.str(idx).zfill(4)+ '_obc.nc'
    ds_tr_east.to_netcdf(fnam,unlimited_dims='time',format='NETCDF3_CLASSIC')

    temp_west = regrid_west_tr(dft['temp'])
    salt_west = regrid_west_tr(dft['salt'])
    temp_west = temp_west.ffill(dim='locations').ffill(dim='st_ocean')
    salt_west = salt_west.ffill(dim='locations').ffill(dim='st_ocean')
    ds_tr_west = xr.Dataset({'temp':temp_west,'salt':salt_west})
    ds_tr_west.time.encoding['calendar']='gregorian'
    fnam=root_folder + 'tracers_west' + np.str(idx).zfill(4)+ '_obc.nc'
    ds_tr_west.to_netcdf(fnam,unlimited_dims='time',format='NETCDF3_CLASSIC')
    # interpolate ssh
    ssh_north = regrid_north_tr(dft['ssh'])
    ds_ssh_north = xr.Dataset({'ssh':ssh_north})
    ds_ssh_north.time.encoding['calendar'] = "gregorian"
    fnam=root_folder + 'ssh_north' + np.str(idx).zfill(4)+ '_obc.nc'
    ds_ssh_north.to_netcdf(fnam,unlimited_dims='time',format='NETCDF3_CLASSIC')
    ssh_south = regrid_south_tr(dft['ssh'])
    ds_ssh_south = xr.Dataset({'ssh':ssh_south})
    ds_ssh_south.time.encoding['calendar'] = "gregorian"
    fnam=root_folder + 'ssh_south' + np.str(idx).zfill(4)+ '_obc.nc'
    ds_ssh_south.to_netcdf(fnam,unlimited_dims='time',format='NETCDF3_CLASSIC')
    ssh_east = regrid_east_tr(dft['ssh'])
    ds_ssh_east = xr.Dataset({'ssh':ssh_east})
    ds_ssh_east.time.encoding['calendar'] = "gregorian"
    fnam=root_folder + 'ssh_east' + np.str(idx).zfill(4)+ '_obc.nc'
    ds_ssh_east.to_netcdf(fnam,unlimited_dims='time',format='NETCDF3_CLASSIC')
    ssh_west = regrid_west_tr(dft['ssh'])
    ds_ssh_west = xr.Dataset({'ssh':ssh_west})
    ds_ssh_west.time.encoding['calendar'] = "gregorian"
    fnam=root_folder + 'ssh_west' + np.str(idx).zfill(4)+ '_obc.nc'
    ds_ssh_west.to_netcdf(fnam,unlimited_dims='time',format='NETCDF3_CLASSIC')
    # interpolate true u and v velocities
    u_north_r = regrid_north_uv(dfuv['u'])
    v_north_r = regrid_north_uv(dfuv['v'])
    u_south_r = regrid_south_uv(dfuv['u'])
    v_south_r = regrid_south_uv(dfuv['v'])
    u_east_r = regrid_east_uv(dfuv['u'])
    v_east_r = regrid_east_uv(dfuv['v'])
    u_west_r = regrid_west_uv(dfuv['u'])
    v_west_r = regrid_west_uv(dfuv['v'])
    # rotate back to model grid
    u_north,v_north=apply_rotation_transpose(u_north_r,v_north_r,ds_regional.angle_dx.isel(nyp=ds_regional.nyp[-1]),time_slice=None)
    u_north = u_north.fillna(0)
    v_north = v_north.fillna(0)
    ds_uv_north = xr.Dataset({'u':u_north,'v':v_north},coords={'lon':north.lon,'lat':north.lat,'st_ocean':dfuv.st_ocean})
    ds_uv_north = ds_uv_north.rename({'i':'locations'})
    ds_uv_north.time.encoding['calendar']='gregorian'
    fnam=root_folder + 'uv_north' + np.str(idx).zfill(4)+ '_obc.nc'
    ds_uv_north.to_netcdf(fnam,unlimited_dims='time',format='NETCDF3_CLASSIC')

    u_south,v_south=apply_rotation_transpose(u_south_r,v_south_r,ds_regional.angle_dx.isel(nyp=ds_regional.nyp[0]),time_slice=None)
    u_south = u_south.fillna(0)
    v_south = v_south.fillna(0)
    ds_uv_south = xr.Dataset({'u':u_south,'v':v_south},coords={'lon':south.lon,'lat':south.lat,'st_ocean':dfuv.st_ocean})
    ds_uv_south = ds_uv_south.rename({'i':'locations'})
    ds_uv_south.time.encoding['calendar']='gregorian'
    fnam=root_folder + 'uv_south' + np.str(idx).zfill(4)+ '_obc.nc'
    ds_uv_south.to_netcdf(fnam,unlimited_dims='time',format='NETCDF3_CLASSIC')

    u_east,v_east=apply_rotation_transpose(u_east_r,v_east_r,ds_regional.angle_dx.isel(nxp=ds_regional.nxp[-1]),time_slice=None)
    u_east = u_east.fillna(0)
    v_east = v_east.fillna(0)
    ds_uv_east = xr.Dataset({'u':u_east,'v':v_east},coords={'lon':east.lon,'lat':east.lat,'st_ocean':dfuv.st_ocean})
    ds_uv_east = ds_uv_east.rename({'i':'locations'})
    ds_uv_east.time.encoding['calendar']='gregorian'
    fnam=root_folder + 'uv_east' + np.str(idx).zfill(4)+ '_obc.nc'
    ds_uv_east.to_netcdf(fnam,unlimited_dims='time',format='NETCDF3_CLASSIC')

    u_west,v_west=apply_rotation_transpose(u_west_r,v_west_r,ds_regional.angle_dx.isel(nxp=ds_regional.nxp[0]),time_slice=None)
    u_west = u_west.fillna(0)
    v_west = v_west.fillna(0)
    ds_uv_west = xr.Dataset({'u':u_west,'v':v_west},coords={'lon':west.lon,'lat':west.lat,'st_ocean':dfuv.st_ocean})
    ds_uv_west = ds_uv_west.rename({'i':'locations'})
    ds_uv_west.time.encoding['calendar']='gregorian'
    fnam=root_folder + 'uv_west' + np.str(idx).zfill(4)+ '_obc.nc'
    ds_uv_west.to_netcdf(fnam,unlimited_dims='time',format='NETCDF3_CLASSIC')





