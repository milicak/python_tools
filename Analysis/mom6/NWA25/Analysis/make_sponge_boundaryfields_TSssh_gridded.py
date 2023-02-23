import os
import numpy as np
import datetime
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
from os.path import exists

root_folder = '/Volumes/A1/workdir/milicak/datasets/MOM6/NWA25/'
mom_dir = '/Volumes/A1/workdir/milicak/datasets/MOM6/NWA25/'
path_regional_grid = mom_dir + './ocean_hgrid.nc'
df2 = xr.open_dataset(path_regional_grid)
lon_rho = np.copy(df2['x'][1::2,1::2])
lat_rho = np.copy(df2['y'][1::2,1::2])
nj,ni = lon_rho.shape
ds2 = df2['x'][1::2,1::2]
ds2 = ds2.to_dataset(name='lon')
ds2['lat']=df2['y'][1::2,1::2]
ds2 = ds2.rename_dims({'nxp': 'x','nyp': 'y'})

encodings = {'time': {'zlib': False, '_FillValue': False},
            'lon': {'zlib': False, '_FillValue': 1e20},
            'lat': {'zlib': False, '_FillValue': 1e20},
            'depth': {'_FillValue': 1e20},
            'temp': {'_FillValue': 1e20},
            'salt': {'_FillValue': 1e20},
            'ssh': {'_FillValue': 1e20},
            }
 # Also fix the time encoding
encodings['time'].update({'dtype':'float64', 'calendar': 'gregorian', 'units': 'hours since 1900-01-01 00:00:00'})

dd1 = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/MOM6/NWA25/OBC_mean/obc_ts_north_1996_mean.nc')
dd2 = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/MOM6/NWA25/OBC_mean/obc_ts_south_1996_mean.nc')
dd3 = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/MOM6/NWA25/OBC_mean/obc_ts_east_1996_mean.nc')

dd1 = dd1.isel(time=slice(0,360,30))
temp1 = dd1.temp_segment_001[:,:,0,1::2]
salt1 = dd1.salt_segment_001[:,:,0,1::2]
dd2 = dd2.isel(time=slice(0,360,30))
temp2 = dd2.temp_segment_002[:,:,0,1::2]
salt2 = dd2.salt_segment_002[:,:,0,1::2]
dd3 = dd3.isel(time=slice(0,360,30))
temp3 = dd3.temp_segment_003[:,:,1::2,0]
salt3 = dd3.salt_segment_003[:,:,1::2,0]

dd1 = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/MOM6/NWA25/OBC_mean/obc_ssh_north_1996_mean.nc')
dd2 = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/MOM6/NWA25/OBC_mean/obc_ssh_south_1996_mean.nc')
dd3 = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/MOM6/NWA25/OBC_mean/obc_ssh_east_1996_mean.nc')

dd1 = dd1.isel(time=slice(0,360,30))
ssh1 = dd1.ssh_segment_001[:,0,1::2]
dd2 = dd2.isel(time=slice(0,360,30))
ssh2 = dd2.ssh_segment_002[:,0,1::2]
dd3 = dd3.isel(time=slice(0,360,30))
ssh3 = dd3.ssh_segment_003[:,1::2,0]

temp = np.zeros((12,50,1690,1550))
salt = np.zeros((12,50,1690,1550))
ssh = np.zeros((12,1690,1550))
# total grid points to sponge
grdpoints = 20
for ind in range(1,grdpoints+2):
    print(ind)
    temp[:,:,:,-ind] = temp3
    temp[:,:,ind,:] = temp2
    temp[:,:,-ind,:] = temp1
    salt[:,:,:,-ind] = salt3
    salt[:,:,ind,:] = salt2
    salt[:,:,-ind,:] = salt1
    ssh[:,:,-ind] = ssh3
    ssh[:,ind,:] = ssh2
    ssh[:,-ind,:] = ssh1


ds_temp = xr.DataArray(temp,coords={'time':temp3.time,'depth':np.copy(temp3.z_l),'lat':np.copy(temp3.lat),'lon':np.copy(temp2.lon)},dims=["time", "depth","lat","lon"])
ds_salt = xr.DataArray(salt,coords={'time':temp3.time,'depth':np.copy(temp3.z_l),'lat':np.copy(temp3.lat),'lon':np.copy(temp2.lon)},dims=["time", "depth","lat","lon"])
ds_ssh = xr.DataArray(ssh,coords={'time':temp3.time,'lat':np.copy(temp3.lat),'lon':np.copy(temp2.lon)},dims=["time", "lat","lon"])

df1 = ds_temp.to_dataset(name='temp')
df2 = ds_salt.to_dataset(name='salt')
df3 = ds_ssh.to_dataset(name='ssh')
# df['salt'] = ds_salt
# df['ssh'] = ds_ssh

out_file = root_folder + 'obc_field_sponge_temp.nc'
df1.to_netcdf(
    out_file,
    unlimited_dims=['time'],
    # encoding=encodings,
    engine='netcdf4')

out_file = root_folder + 'obc_field_sponge_salt.nc'
df2.to_netcdf(
    out_file,
    unlimited_dims=['time'],
    # encoding=encodings,
    engine='netcdf4')

out_file = root_folder + 'obc_field_sponge_ssh.nc'
df3.to_netcdf(
    out_file,
    unlimited_dims=['time'],
    # encoding=encodings,
    engine='netcdf4')
# ncatted -a cartesian_axis,depth,c,c,"Z" SODA_TSssh_sponge_v2.nc


