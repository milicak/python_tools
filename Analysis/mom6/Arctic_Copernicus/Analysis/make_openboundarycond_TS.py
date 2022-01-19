import os
import numpy as np
import xesmf as xe
import xarray as xr
import scipy.io
from scipy.io import savemat
from scipy.io import loadmat
from mpl_toolkits.basemap import Basemap, shiftgrid
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


#Note that parent grid uv values are symmetric
# path_parent_grid='/net2/mjh/ipynb/OM4_025/c192_mosaic/ocean_hgrid.nc'
# parent_grid=open_grid(path_parent_grid)
path_regional_grid='./ocean_hgrid.nc'
regional_grid=open_grid(path_regional_grid)

dsr_topo=xr.open_dataset('ocean_topog.nc')
dsr_topo = xr.merge([dsr_topo, regional_grid['h']])
dsr_topo.depth.plot(vmax=-50.,vmin=250.,cmap='gist_gray')
txt=plt.title('Regional Domain')


variable = 'templvl'

ds = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic/NOIIAJRAOC20TR_TL319_tn14_20190710.micom.hm.1958-01.nc')
gr1 = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic/NorESM_grid.nc')
mask = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/ocean_mask.nc')

df = flood_kara(ds[variable], xdim='x', ydim='y', zdim='depth')
# df = df[0,0,:,:]
df = df.to_dataset(name='temp')
df['lon'] = gr1['plon']
df['lat'] = gr1['plat']

# output grid info
df2 = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/ocean_hgrid.nc')
lon_rho = np.copy(df2['x'][1::2,1::2])
lat_rho = np.copy(df2['y'][1::2,1::2])
nj,ni = lon_rho.shape
ds2 = df2['x'][1::2,1::2]
ds2 = ds2.to_dataset(name='lon')
ds2['lat']=df2['y'][1::2,1::2]
ds2 = ds2.rename_dims({'nxp': 'x','nyp': 'y'})

# for regridding
temp = np.copy(df.temp)
dstemp = xr.DataArray(temp, coords=[df.time, df.depth, df.y, df.x],
                            dims=["time","depth","y","x"])
df['temp'] = dstemp

# build regridder
regridder = xe.Regridder(df, ds2, 'nearest_s2d', reuse_weights=True)

