import os
import numpy as np
import xesmf as xe
import xarray as xr
import scipy.io
from scipy.io import savemat
from scipy.io import loadmat
import matplotlib.colors as colors
from scipy.signal import medfilt2d
import netCDF4
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.path import Path
#for interpolation
from scipy.spatial import cKDTree
from HCtFlood.kara import flood_kara

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
    # deg_rad=np.pi/180.
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

    return u,v

root_folder = '/okyanus/users/milicak/dataset/SODA3_12_21/'
fname1 = 'SODA_monthly_sponge_1980.nc'
ds = xr.open_dataset(root_folder + fname1)
mom_dir = '/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/'
path_regional_grid = mom_dir + './ocean_hgrid.nc'
# Define boundaries for regional domain model
regional_grid = open_grid(path_regional_grid)
ds_regional = regional_grid['ds']
dsangle = ds_regional.angle_dx[:-1:2,:-1:2]


variable = 'u'
dft = ds[variable][:,:,200:,:]
dft = dft.to_dataset(name=variable)
dft = dft.rename({'xu_ocean': 'lon', 'yu_ocean': 'lat'})
# dft = dft.fillna(0)
dft = dft.rename({'st_ocean':'depth'})
dft = dft.ffill(dim='depth')
dft = flood_kara(dft[variable], xdim='lon', ydim='lat', zdim='depth')
dft.load()
dft = dft.to_dataset(name=variable)
# salinity
variable = 'v'
dfs = ds[variable][:,:,200:,:]
dfs = dfs.to_dataset(name=variable)
dfs = dfs.rename({'xu_ocean': 'lon', 'yu_ocean': 'lat'})
# dfs = dfs.fillna(0)
dfs = dfs.rename({'st_ocean':'depth'})
dfs = dfs.ffill(dim='depth')
dfs = flood_kara(dfs[variable], xdim='lon', ydim='lat', zdim='depth')
dfs.load();
dfs = dfs.to_dataset(name=variable)

df2 = xr.open_dataset(path_regional_grid)
lon_rho = np.copy(df2['x'][1::2,1::2])
lat_rho = np.copy(df2['y'][1::2,1::2])
nj,ni = lon_rho.shape
ds2 = df2['x'][1::2,1::2]
ds2 = ds2.to_dataset(name='lon')
ds2['lat']=df2['y'][1::2,1::2]
ds2 = ds2.rename_dims({'nxp': 'x','nyp': 'y'})

# build regridder
# regridder = xe.Regridder(df, ds2, 'nearest_s2d',reuse_weights=True)

#apply regridder
# dr_out1 = regridder(df[variable])

#apply regridder
# dr_out2 = regridder(df[variable])

df = dfs.drop_vars('v')
df = df.drop('depth')

# build regridder
regridder2 = xe.Regridder(df, ds2, 'bilinear',reuse_weights=True)

variable = 'ssh'
dfssh = ds[variable][:,200:,:]
dfssh = dfssh.to_dataset(name=variable)
dfssh = dfssh.rename({'xt_ocean': 'lon', 'yt_ocean': 'lat'})
dfssh = dfssh.fillna(0)
dr_out3 = regridder2(dfssh[variable])
dr_out3 = dr_out3.where(dr_out3!=0)
dr_out3 = dr_out3.to_dataset(name=variable)
dr_out3 = flood_kara(dr_out3[variable], xdim='x', ydim='y')
dr_out3.load();
dr_out3 = dr_out3.isel(z=0);
dr_out3 = dr_out3.drop('z');

uvel = np.zeros((13,50,1300,2000))
vvel = np.zeros((13,50,1300,2000))
for kind in range(0,50):
    print(kind)
    tmp = dft.isel(depth=kind)
    tmp = tmp.drop('depth')
    smp = dfs.isel(depth=kind)
    smp = smp.drop('depth')
    dr_out1 = regridder2(tmp)
    dr_out1 = dr_out1.where(dr_out1!=0)
    dr_out1 = dr_out1.interpolate_na(dim="y", method="nearest",
                            fill_value="extrapolate")
    dr_out2 = regridder2(smp)
    dr_out2 = dr_out2.where(dr_out2!=0)
    dr_out2 = dr_out2.interpolate_na(dim="y", method="nearest",
                            fill_value="extrapolate")
    u_north,v_north = apply_rotation_transpose(dr_out1.u,dr_out2.v,dsangle,time_slice=None)
    uvel[:,kind,:,:] = np.copy(u_north)
    vvel[:,kind,:,:] = np.copy(v_north)



# make a dataarray
ssh = np.copy(dr_out3)

data_vars = {'uvel':(['time','depth','y','x'], uvel,),
             'vvel':(['time','depth','y','x'], vvel,),
             'lon':(['y','x'], np.copy(ds2.lon),),
             'lat':(['y','x'], np.copy(ds2.lat),),
             'depth':(['depth'], np.copy(dft.depth),),
            }


# define coordinates
coords = {'time': (['time'], dr_out3.time),
          # 'depth': (['depth'], dft.depth),
          'y': (['y'], dr_out3.y),
          'x': (['x'], dr_out3.x),
         }

# create dataset
dr_out = xr.Dataset(data_vars=data_vars,
                coords=coords)


encodings = {'time': {'zlib': False, '_FillValue': False},
            'lon': {'zlib': False, '_FillValue': 1e20},
            'lat': {'zlib': False, '_FillValue': 1e20},
            'depth': {'_FillValue': 1e20},
            'uvel': {'_FillValue': 1e20},
            'vvel': {'_FillValue': 1e20},
            }

# Make sure time has the right units and datatype
# otherwise it will become an int and MOM will fail.
encodings['time'].update({
        'units': 'days since 1950-01-01',
        'dtype': np.float,
        'calendar': 'NOLEAP'})

encodings['depth'].update({
        # 'long_name': 'tcell zstar depth',
        # 'units': 'meters',
        # 'positive': 'up',
        # 'axis': 'Z',
        # 'cartesian_axis': 'Z'
        })

out_file = mom_dir + 'SODA_UV_sponge_v2.nc'

dr_out.to_netcdf(
    out_file,
    unlimited_dims=['time'],
    # format='NETCDF3_64BIT',
    encoding=encodings,
    engine='netcdf4')

dr_out.close()
# ncatted -a cartesian_axis,depth,c,c,"Z" SODA_TSssh_sponge_v2.nc
# ncatted -a axis,depth,c,c,"Z" SODA_TSssh_sponge_v2.nc
# ncatted -a positive,depth,c,c,"up" SODA_TSssh_sponge_v2.nc
# ncatted -a units,depth,c,c,"meters" SODA_TSssh_sponge_v2.nc
# ncks -C -O -x -v x,y SODA_TSssh_sponge_v2.nc
