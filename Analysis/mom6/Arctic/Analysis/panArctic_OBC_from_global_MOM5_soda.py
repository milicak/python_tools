import xarray as xr
import xesmf as xe
from xgcm import Grid
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


def velocity_at_corners(ds_u,ds_v):
    nxp=ds_u.nxp[-1].data+1;nyp=ds_v.nyp[-1].data+1
    #upper-right q points
    u_q=0.5*(ds_u.uo+ds_u.uo.roll(roll_coords='yh',yh=-1)).isel(xq=slice(1,nxp))
    #upper-right q points
    v_q=0.5*(ds_v.vo+ds_v.vo.roll(roll_coords='xh',xh=-1)).isel(yq=slice(1,nyp))
    ds_uvq = xr.Dataset({'uo':u_q,'vo':v_q},coords={'time':ds_u.time,'lon':parent_grid['q'].x,'lat':parent_grid['q'].y,'angle_dx':parent_grid['q'].angle_dx})
    return ds_uvq


def velocity_at_corners_nonsym(ds_u,ds_v):
    nxp=ds_u.nxp[-1].data+1;nyp=ds_v.nyp[-1].data+1
    #upper-right q points
    u_q=0.5*(ds_u.uo+ds_u.uo.roll(roll_coords='yh',yh=-1)).isel(xq=slice(0,nxp))
    #upper-right q points
    v_q=0.5*(ds_v.vo+ds_v.vo.roll(roll_coords='xh',xh=-1)).isel(yq=slice(0,nyp))
    ds_uvq = xr.Dataset({'uo':u_q,'vo':v_q},coords={'time':ds_u.time,'lon':parent_grid['q'].x,'lat':parent_grid['q'].y,'angle_dx':parent_grid['q'].angle_dx})
    return ds_uvq



def velocity_at_corners_grid(ds_u,ds_v,grid):
    nxp=ds_u.nxp[-1].data+1;nyp=ds_v.nyp[-1].data+1
    #upper-right q points
    u_q = grid.interp(ds_u.uo, 'Y', boundary='fill')
    #upper-right q points
    v_q = grid.interp(ds_v.vo, 'X', boundary='fill')
    ds_uvq = xr.Dataset({'uo':u_q,'vo':v_q},coords={'time':ds_u.time,'lon':parent_grid['q'].x,'lat':parent_grid['q'].y,'angle_dx':parent_grid['q'].angle_dx})
    return ds_uvq



path_parent_grid='/archive/milicak/dataset/MOM6/OM4_05/mosaic_ocean.v20180227.unpacked/ocean_hgrid.nc'
parent_grid=open_grid(path_parent_grid)
path_regional_grid='/archive/milicak/dataset/MOM6/Arctic/ocean_hgrid.nc'
regional_grid=open_grid(path_regional_grid)
regional_grid['ds']=regional_grid['ds'].load()


# plot the angles
regional_grid['ds'].angle_dx.plot()

# load the topo file
dsr_topo=xr.open_dataset('/archive/milicak/dataset/MOM6/Arctic/topog.nc')
dsr_topo = xr.merge([dsr_topo, regional_grid['h']])
#dsr_topo.depth.plot(vmax=-50.,vmin=250.,cmap='gist_gray')
plt.pcolormesh(dsr_topo.depth.data,vmin=-50,vmax=250,cmap=plt.cm.gist_gray)
txt=plt.title('Regional Domain')

path_model_data='/archive/milicak/MOM6-examples/Projects/patchy_NA/work_ctrl/OUT/19350101.ocean_month_z.nc'
df = xr.open_dataset(path_model_data)
grid = Grid(df, coords={'X': {'center': 'xh', 'right': 'xq'},
                        'Y': {'center': 'yh', 'right': 'yq'},
                        'Z': {'inner': 'z_l', 'outer': 'z_i'} }, periodic=['X'])
fields=[{'temp':'h'},{'salt':'h'},{'ssh':'h'},{'u':'Cu'},{'v':'Cv'}]
# fields=[{'thetao':'h'},{'so':'h'},{'ssh':'h'},{'uo':'Cu'},{'vo':'Cv'}]
# fields=[{'thetao':'h'}]
model_data = open_dataset(path_model_data,fields,parent_grid)

# u,v variables
ds_u=model_data['ds_u'];ds_v=model_data['ds_v']
# model_data['ds_uv']=velocity_at_corners(ds_u,ds_v)
model_data['ds_uv'] = velocity_at_corners_grid(ds_u,ds_v,grid)


ds_uvq = model_data['ds_uv']
uq_rot,vq_rot = apply_rotation(ds_uvq.uo,ds_uvq.vo,ds_uvq.angle_dx,time_slice=None)
ds_uvq_r = xr.Dataset({'u':uq_rot,'v':vq_rot},coords={'time':ds_uvq.time,'lon':ds_uvq.lon,'lat':ds_uvq.lat})
model_data['ds_uv_r'] = ds_uvq_r

# Define open boundary conditions
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


#Load tracers
ds_tr=model_data['ds_tr']
# ds_tr = xr.Dataset({'temp':ds_tr.temp,'salt':ds_tr.salt},coords={'lon':ds_tr.x,'lat':ds_tr.y})
ds_tr = xr.Dataset({'temp':ds_tr.thetao,'salt':ds_tr.so},coords={'lon':ds_tr.x,'lat':ds_tr.y})

# Compute weights
# Using nearest neighbor - other options could be used here , e.g. bilinear.
regrid_north_uv = xe.Regridder(ds_uvq_r, north, 'nearest_s2d',
                               locstream_out=True, periodic=False, filename='regrid_north_uv.nc',reuse_weights=False)
regrid_south_uv = xe.Regridder(ds_uvq_r, south, 'nearest_s2d',
                               locstream_out=True, periodic=False, filename='regrid_south_uv.nc',reuse_weights=False)
regrid_east_uv = xe.Regridder(ds_uvq_r, east, 'nearest_s2d',
                               locstream_out=True, periodic=False, filename='regrid_east_uv.nc',reuse_weights=False)
regrid_west_uv = xe.Regridder(ds_uvq_r, west, 'nearest_s2d',
                               locstream_out=True, periodic=False, filename='regrid_west_uv.nc',reuse_weights=False)

regrid_north_tr = xe.Regridder(ds_tr, north, 'nearest_s2d',
                               locstream_out=True, periodic=False, filename='regrid_north_tr.nc',reuse_weights=False)
regrid_south_tr = xe.Regridder(ds_tr, south, 'nearest_s2d',
                               locstream_out=True, periodic=False, filename='regrid_south_tr.nc',reuse_weights=False)
regrid_east_tr = xe.Regridder(ds_tr, east, 'nearest_s2d',
                               locstream_out=True, periodic=False, filename='regrid_east_tr.nc',reuse_weights=False)
regrid_west_tr = xe.Regridder(ds_tr, west, 'nearest_s2d',
                               locstream_out=True, periodic=False, filename='regrid_west_tr.nc',reuse_weights=False)



# compute u,v boundaries
u_north_r = regrid_north_uv(ds_uvq_r['u'])
v_north_r = regrid_north_uv(ds_uvq_r['v'])
u_south_r = regrid_south_uv(ds_uvq_r['u'])
v_south_r = regrid_south_uv(ds_uvq_r['v'])
u_west_r = regrid_west_uv(ds_uvq_r['u'])
v_west_r = regrid_west_uv(ds_uvq_r['v'])
u_east_r = regrid_east_uv(ds_uvq_r['u'])
v_east_r = regrid_east_uv(ds_uvq_r['v'])

u_north,v_north=apply_rotation_transpose(u_north_r,v_north_r,ds_regional.angle_dx.isel(nyp=ds_regional.nyp[-1]),time_slice=None)
ds_uv_north = xr.Dataset({'u':u_north,'v':v_north},coords={'lon':north.lon,'lat':north.lat,'z_l':ds_uvq.z_l})
ds_uv_north.time.attrs['calendar']='gregorian'
fnam='uv_north.nc'
ds_uv_north.to_netcdf(fnam,unlimited_dims='time',format='NETCDF3_CLASSIC')
u_south,v_south=apply_rotation_transpose(u_south_r,v_south_r,ds_regional.angle_dx.isel(nyp=ds_regional.nyp[0]),time_slice=None)
ds_uv_south = xr.Dataset({'u':u_south,'v':v_south},coords={'lon':south.lon,'lat':south.lat,'z_l':ds_uvq.z_l})
fnam='uv_south.nc'
ds_uv_south.time.attrs['calendar']='gregorian'
ds_uv_south.to_netcdf(fnam,unlimited_dims='time',format='NETCDF3_CLASSIC')
u_west,v_west=apply_rotation_transpose(u_west_r,v_west_r,ds_regional.angle_dx.isel(nxp=ds_regional.nxp[0]),time_slice=None)
ds_uv_west = xr.Dataset({'u':u_west,'v':v_west},coords={'lon':west.lon,'lat':west.lat,'z_l':ds_uvq.z_l})
fnam='uv_west.nc'
ds_uv_west.time.attrs['calendar']='gregorian'
ds_uv_west.to_netcdf(fnam,unlimited_dims='time',format='NETCDF3_CLASSIC')
u_east,v_east=apply_rotation_transpose(u_east_r,v_east_r,ds_regional.angle_dx.isel(nxp=ds_regional.nxp[-1]),time_slice=None)
ds_uv_east = xr.Dataset({'u':u_east,'v':v_east},coords={'lon':east.lon,'lat':east.lat,'z_l':ds_uvq.z_l})
fnam='uv_east.nc'
ds_uv_east.time.attrs['calendar']='gregorian'
ds_uv_east.to_netcdf(fnam,unlimited_dims='time',format='NETCDF3_CLASSIC')


# compute temp/salt boundaries
temp_north = regrid_north_tr(ds_tr['temp'])
salt_north = regrid_north_tr(ds_tr['salt'])
ds_tr_north = xr.Dataset({'temp':temp_north,'salt':salt_north})
ds_tr_north.time.attrs['calendar']='gregorian'
fnam='tracers_north.nc'
ds_tr_north.to_netcdf(fnam,unlimited_dims='time',format='NETCDF3_CLASSIC')
temp_south = regrid_south_tr(ds_tr['temp'])
salt_south = regrid_south_tr(ds_tr['salt'])
ds_tr_south = xr.Dataset({'temp':temp_south,'salt':salt_south})
ds_tr_south.time.attrs['calendar']='gregorian'
fnam='tracers_south.nc'
ds_tr_south.to_netcdf(fnam,unlimited_dims='time',format='NETCDF3_CLASSIC')
temp_west = regrid_west_tr(ds_tr['temp'])
salt_west = regrid_west_tr(ds_tr['salt'])
ds_tr_west = xr.Dataset({'temp':temp_west,'salt':salt_west})
ds_tr_west.time.attrs['calendar']='gregorian'
fnam='tracers_west.nc'
ds_tr_west.to_netcdf(fnam,unlimited_dims='time',format='NETCDF3_CLASSIC')
temp_east = regrid_east_tr(ds_tr['temp'])
salt_east = regrid_east_tr(ds_tr['salt'])
ds_tr_east = xr.Dataset({'temp':temp_east,'salt':salt_east})
ds_tr_east.time.attrs['calendar']='gregorian'
fnam='tracers_east.nc'
ds_tr_east.to_netcdf(fnam,unlimited_dims='time',format='NETCDF3_CLASSIC')


