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

mom_dir = '/okyanus/users/milicak/dataset/MOM6/NA12/'
path_regional_grid = mom_dir + './ocean_hgrid.nc'
df2 = xr.open_dataset(path_regional_grid)
lon_rho = np.copy(df2['x'][1::2,1::2])
lat_rho = np.copy(df2['y'][1::2,1::2])
nj,ni = lon_rho.shape
ds2 = df2['x'][1::2,1::2]
ds2 = ds2.to_dataset(name='lon')
ds2['lat']=df2['y'][1::2,1::2]
ds2 = ds2.rename_dims({'nxp': 'x','nyp': 'y'})

root_folder = '/okyanus/users/milicak/dataset/SODA3_12_21/monthlyfiles/'
f = 'soda3.12.2_mn_ocean_reg'

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

for year in range(2000,2018):
    print(year)
    if year==2017:
        current = xr.open_dataset(f"{root_folder}/{f}_{year}.nc")
        previous = xr.open_dataset(f"{root_folder}/{f}_{year-1}.nc").isel(time=-1)
        out = xr.concat([previous, current], dim="time")
    else:
        current = xr.open_dataset(f"{root_folder}/{f}_{year}.nc")
        next_data = xr.open_dataset(f"{root_folder}/{f}_{year+1}.nc").isel(time=0)
        previous = xr.open_dataset(f"{root_folder}/{f}_{year-1}.nc").isel(time=-1)
        out = xr.concat([previous, current, next_data], dim="time")


    out = out.rename({'xt_ocean': 'lon', 'yt_ocean': 'lat','st_ocean':'depth'})
    # swap dimensions
    out2 = out.temp.transpose("time","depth","lat","lon")
    out2 = out2.to_dataset(name='temp')
    out2['salt'] = out.salt.transpose("time","depth","lat","lon")
    out2['ssh'] = out.ssh.transpose("time","lat","lon")

    if year==2017:
        newtime = xr.Dataset({'time': datetime.datetime(2018, 1, 16)})
        out3 = out2.isel(time=-1)
        out3['time'] = newtime.time
        out2 = xr.concat([out2, out3], dim="time")


    # build regridder
    if exists('patch_330x720_1844x1678_peri.nc'):
        regridder2 = xe.Regridder(out2, ds2, 'patch',reuse_weights=True, periodic=True)
    else:
        regridder2 = xe.Regridder(out2, ds2, 'patch',reuse_weights=False, periodic=True)


    dft = flood_kara(out2['temp'], xdim='lon', ydim='lat', zdim='depth')
    dft.load();
    dft = dft.to_dataset(name='temp')
    dst = flood_kara(out2['salt'], xdim='lon', ydim='lat', zdim='depth')
    dst.load();
    dft['salt'] = dst
    dzt = flood_kara(out2['ssh'], xdim='lon', ydim='lat')
    dzt.load();
    dzt = dzt.to_dataset(name='ssh')
    dzt = dzt.ssh.isel(z=0)
    dzt = dzt.drop('z')
    dft['ssh'] = dzt

    dr_out = regridder2(dft)
    out_file = '/okyanus/users/milicak/dataset/MOM6/NA12/SODA_sponge_' + str(year) + '.nc'
    dr_out.to_netcdf(
        out_file,
        unlimited_dims=['time'],
        # format='NETCDF3_64BIT',
        encoding=encodings,
        engine='netcdf4')



# ncatted -a cartesian_axis,depth,c,c,"Z" SODA_TSssh_sponge_v2.nc
