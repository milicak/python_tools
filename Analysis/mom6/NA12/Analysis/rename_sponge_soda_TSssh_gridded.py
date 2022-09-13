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

df = xr.open_dataset('~/dataset/MOM6/NA12/SODA_sponge_1996.nc')

encodings = {'time': {'zlib': False, '_FillValue': False},
            'lon': {'zlib': False, '_FillValue': 1e20},
            'lat': {'zlib': False, '_FillValue': 1e20},
            'z_l': {'_FillValue': 1e20},
            'z_i': {'_FillValue': 1e20},
            'temp': {'_FillValue': 1e20},
            'salt': {'_FillValue': 1e20},
            'ssh': {'_FillValue': 1e20},
            }
 # Also fix the time encoding
encodings['time'].update({'dtype':'float64', 'calendar': 'gregorian', 'units': 'hours since 1900-01-01 00:00:00'})
df = df.rename({'x': 'xh', 'y': 'yh','depth':'z_l'})
df.z_l.attrs['positive'] = 'down'
dz = np.copy(df.z_l[1:])-np.copy(df.z_l[:-1])
zw = np.concatenate([np.array([0],float),dz]).cumsum()
zw = np.concatenate([zw,np.array([8500],float)])
df['z_i'] = zw
year = 1996
out_file = '/okyanus/users/milicak/dataset/MOM6/NA12/SODA_sponge_' + str(year) + '_new.nc'
df.to_netcdf(
    out_file,
    unlimited_dims=['time'],
    # format='NETCDF3_64BIT',
    encoding=encodings,
    engine='netcdf4')



# ncatted -a cartesian_axis,depth,c,c,"Z" SODA_TSssh_sponge_v2.nc
