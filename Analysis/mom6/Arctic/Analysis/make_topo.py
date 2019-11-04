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


df = xr.open_dataset('roms_grid_orig.nc')
lon_roms = np.copy(df.lon_rho)
lat_roms = np.copy(df.lat_rho)

nj,ni = df.lon_rho.shape

# Read and interpolate TOPO file into the grid
gr = xr.open_dataset('/okyanus/users/milicak/dataset/world_bathy/GEBCO_2014_1D.nc')
lon_gebco = np.linspace(gr.x_range[0],gr.x_range[1],gr.dimension[0],endpoint=True);
lat_gebco = np.linspace(gr.y_range[0],gr.y_range[1],gr.dimension[1],endpoint=True);
depth_gebco = np.reshape(np.copy(gr.z), (np.copy(gr.dimension[1]), np.copy(gr.dimension[0])))
lat_gebco = lat_gebco[15000:]
depth_gebco = np.flipud(depth_gebco)
depth_gebco = depth_gebco[15000:,:]

# Create xarray dataset for depth_gebco
foo = xr.Dataset({'depth_g':(['lat','lon'],  depth_gebco)})
foo = foo.assign_coords(lat=lat_gebco,lon=lon_gebco)

dd = foo.interp(lat=df.lat_rho,lon=df.lon_rho)

aa = np.copy(dd.depth_g.where(dd.depth_g<0,0))
aa[(aa>-5) & (aa<0)]=-5
aa[aa<-6000]=-6000
hraw = -aa

# for Arctic setup

# Create a topography file
rg = scipy.io.netcdf_file('ocean_topog.nc','w')
# Dimensions
rg.createDimension('nx',ni)
rg.createDimension('ny',nj)
# Variables
hdepth = rg.createVariable('depth','float32',('ny','nx',))
hdepth.units = 'm'
# Values
hdepth[:] = hraw #[0,1:-1,1:-1]
rg.close()

