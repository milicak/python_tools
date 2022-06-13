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


root_folder = '/okyanus/users/milicak/dataset/SODA3_12_21/'
fname1 = 'SODA_monthly_sponge_1980.nc'
ds = xr.open_dataset(root_folder + fname1)
mom_dir = '/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/'
path_regional_grid = mom_dir + './ocean_hgrid.nc'


variable = 'temp'
df = ds[variable][:,:,200:,:]
df = df.to_dataset(name=variable)
df = df.rename({'xt_ocean': 'lon', 'yt_ocean': 'lat'})
df = df.fillna(0)
df = df.rename({'st_ocean':'depth'})
df = df.ffill(dim='depth')

df2 = xr.open_dataset(path_regional_grid)
lon_rho = np.copy(df2['x'][1::2,1::2])
lat_rho = np.copy(df2['y'][1::2,1::2])
nj,ni = lon_rho.shape
ds2 = df2['x'][1::2,1::2]
ds2 = ds2.to_dataset(name='lon')
ds2['lat']=df2['y'][1::2,1::2]
ds2 = ds2.rename_dims({'nxp': 'x','nyp': 'y'})

# build regridder
regridder = xe.Regridder(df, ds2, 'nearest_s2d',reuse_weights=True)

#apply regridder
dr_out1 = regridder(df[variable])

# salinity
variable = 'salt'
df = ds[variable][:,:,200:,:]
df = df.to_dataset(name=variable)
df = df.rename({'xt_ocean': 'lon', 'yt_ocean': 'lat'})
df = df.fillna(0)
df = df.rename({'st_ocean':'depth'})
df = df.ffill(dim='depth')
#apply regridder
dr_out2 = regridder(df[variable])

variable = 'ssh'
df = ds[variable][:,:,:]
df = df.to_dataset(name=variable)
df = df.rename({'xt_ocean': 'lon', 'yt_ocean': 'lat'})
df = df.fillna(0)

# build regridder
# regridder2 = xe.Regridder(df, ds2, 'nearest_s2d',reuse_weights=True)
regridder2 = xe.Regridder(df, ds2, 'bilinear',reuse_weights=True)
#apply regridder
dr_out3 = regridder2(df[variable])

encodings = {'time': {'zlib': False, '_FillValue': False},
            'lon': {'zlib': False, '_FillValue': 1e20},
            'lat': {'zlib': False, '_FillValue': 1e20},
            'temp': {'_FillValue': 1e20},
            'salt': {'_FillValue': 1e20},
            'ssh': {'_FillValue': 1e20},
            }

# Make sure time has the right units and datatype
# otherwise it will become an int and MOM will fail.
encodings['time'].update({
        'units': 'days since 1950-01-01',
        'dtype': np.float,
        'calendar': 'NOLEAP'})

dr_out = dr_out1.to_dataset(name='temp')
dr_out['salt'] = dr_out2
dr_out['ssh'] = dr_out3
out_file = mom_dir + 'SODA_TSssh_sponge.nc'

dr_out.to_netcdf(
    out_file,
    unlimited_dims=['time'],
    # format='NETCDF3_64BIT',
    encoding=encodings,
    engine='netcdf4')

dr_out.close()

