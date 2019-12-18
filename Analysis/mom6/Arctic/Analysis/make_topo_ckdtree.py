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


df = xr.open_dataset('roms_grid_orig.nc')
lon_roms = np.copy(df.lon_rho)
lat_roms = np.copy(df.lat_rho)

nj,ni = df.lon_rho.shape

# Read and interpolate TOPO file into the grid
gr = xr.open_dataset('/archive/milicak/dataset/world_bathy/GEBCO_2014_1D.nc')
lon_gebco = np.linspace(gr.x_range[0],gr.x_range[1],gr.dimension[0],endpoint=True);
lat_gebco = np.linspace(gr.y_range[0],gr.y_range[1],gr.dimension[1],endpoint=True);
depth_gebco = np.reshape(np.copy(gr.z), (np.copy(gr.dimension[1]), np.copy(gr.dimension[0])))

# for coarse resolution
# lon_gebco=lon_gebco[::10]
# lat_gebco=lat_gebco[::10]
# depth_gebco=depth_gebco[::10,::10]

depth_gebco = np.flipud(depth_gebco)
lon_gebco,lat_gebco = np.meshgrid(lon_gebco,lat_gebco)

# for Arctic setup
depth_gebco = depth_gebco[15000:,:]
lon_gebco = lon_gebco[15000:,:]
lat_gebco = lat_gebco[15000:,:]

# points = np.column_stack((lon_gebco.flatten(), lat_gebco.flatten()))
# hraw = griddata((lon_gebco.flatten(),lat_gebco.flatten(),depth_gebco.flatten(),(np.copy(df.lon_rho),np.copy(df.lat_rho)))
# hraw = griddata(points,depth_gebco.flatten(),(np.copy(df.lon_rho),np.copy(df.lat_rho)))
# hraw.shape = df.lon_rho.shape

# Use cKDTree
xs, ys, zs = lon_lat_to_cartesian(lon_gebco.flatten(), lat_gebco.flatten())
del zs
xt, yt, zt = lon_lat_to_cartesian(lon_roms.flatten(), lat_roms.flatten())
del zt

# aa = np.transpose(np.vstack((lon_gebco.flatten(),lat_gebco.flatten())))
aa = np.column_stack((lon_gebco.flatten(),lat_gebco.flatten()))
bb = np.column_stack((xs,ys))
tree = cKDTree(aa)
d, inds = tree.query(np.transpose(np.vstack((np.copy(df.lon_rho).flatten(),
                                             np.copy(df.lat_rho).flatten()))),
                     k = 4)
w = 1.0 / d**2
hraw = np.sum(w * depth_gebco.flatten()[inds], axis=1) / np.sum(w, axis=1)
hraw.shape = df.lon_rho.shape

# x, y = np.meshgrid(np.arange(300), np.arange(300)) # make a canvas with coordinates
# x, y = x.flatten(), y.flatten()
# points = np.vstack((x,y)).T
# p = Path(tupVerts) # make a polygon
# grid = p.contains_points(points)
# mask = grid.reshape(300,300) # now you have a mask with points inside a polygon

# Create a topography file
rg = scipy.io.netcdf_file('ocean_topog.nc','w')
# Dimensions
rg.createDimension('ni',ni)
rg.createDimension('nj',nj)
# Variables
hdepth = rg.createVariable('depth','float32',('nj','ni',))
hdepth.units = 'm'
# Values
hdepth[:] = hraw #[0,1:-1,1:-1]
rg.close()

