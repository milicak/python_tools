# This file is designed to be cut and pasted into an ipython --pylab
# session. Otherwise, you'll need to "import np as np" then
# convert "array" to "np.array".
import os
import numpy as np
import xarray as xr
import scipy.io
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.colors as colors
from scipy.signal import medfilt2d
import netCDF4
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.path import Path
#for interpolation
from scipy.spatial import cKDTree


def lon_lat_to_cartesian(lon, lat, R = 1):
    """
    calculates lon, lat coordinates of a point on a sphere with
    radius R
    """
    lon_r = np.radians(lon)
    lat_r = np.radians(lat)
    x =  R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)
    return x,y,z



df = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic/ocean_hgrid.nc')
lon_rho = np.copy(df['x'][1::2,1::2])
lat_rho = np.copy(df['y'][1::2,1::2])
nj,ni = lon_rho.shape

# input data
gr = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/OM4_025/ocean_static.nc')
lon_input = np.copy(gr.geolon)
lat_input = np.copy(gr.geolat)
di = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/OM4_025/INPUT/tidal_amplitude.v20140616.nc');
tide_in = np.copy(di.tideamp)

xs, ys, zs = lon_lat_to_cartesian(lon_input.flatten(), lat_input.flatten())
xt, yt, zt = lon_lat_to_cartesian(lon_rho.flatten(), lat_rho.flatten())

tree = cKDTree(list(zip(xs, ys, zs)))
d, inds = tree.query(list(zip(xt, yt, zt)), k = 4)
w = 1.0 / d**2
tide_out = np.sum(w * tide_in.flatten()[inds], axis=1) / np.sum(w, axis=1)
tide_out.shape = lon_rho.shape

# Create a topography file
rg = scipy.io.netcdf_file('tidal_amplitude.nc','w')
# Dimensions
rg.createDimension('nx',ni)
rg.createDimension('ny',nj)
# Variables

nx = rg.createVariable('nx','float32',('nx',))
nx.units = 'degrees_east'
nx[:] = lon_rho[0,:]
ny = rg.createVariable('ny','float32',('ny',))
ny.units = 'degrees_north'
ny[:] = lat_rho[:,0]

tideamp = rg.createVariable('tideamp','float32',('ny','nx',))
tideamp.units = 'm s-1'
# Values
tideamp[:] = tide_out #[0,1:-1,1:-1]

rg.close()

