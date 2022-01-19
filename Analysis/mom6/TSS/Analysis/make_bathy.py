# This file is designed to be cut and pasted into an ipython --pylab
# session. Otherwise, you'll need to "import np as np" then
# convert "array" to "np.array".
import os
import numpy as np
import scipy.io
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.colors as colors
from scipy.signal import medfilt2d
import netCDF4
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import pyroms
import pyroms_toolbox
from bathy_smoother import *
from matplotlib.path import Path
#for interpolation
from scipy.spatial import cKDTree



hgrd = scipy.io.loadmat('MITgcm_TSS_grid.mat')
lon_rho = np.transpose(hgrd['lon_rho'])
lat_rho = np.transpose(hgrd['lat_rho'])
depth = np.transpose(hgrd['bathy'])
depth = depth[0:-1,0:-1]

nj,ni = lon_rho.shape
nj -=2; ni -=2
print('nj=%i, nj=%i'%(nj,ni))

# Declare shapes
bathy = np.zeros((nj,ni))

# Copy in data from ROMS file
bathy = depth[1:-1,1:-1] # Cell centers (drop outside row and column)

# Create a mosaic file
rg = scipy.io.netcdf_file('/okyanus/users/milicak/dataset/MOM6/TSS/ocean_topog.nc','w')
# Dimensions
rg.createDimension('nx',ni)
rg.createDimension('ny',nj)
rg.createDimension('ntiles',1)
# Variables
hx = rg.createVariable('depth','float32',('ny','nx',))
hx.units = 'm'
# Values
hx[:,:] = bathy
rg.close()



