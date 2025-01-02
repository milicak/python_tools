# This file is designed to be cut and pasted into an ipython --pylab
# session. Otherwise, you'll need to "import np as np" then
# convert "array" to "np.array".
import os
import numpy as np
import xarray as xr
import scipy.io
import netCDF4
# first run the following

# cd ~/python_libs/ocean_model_topog_generator
# OMtopogen/create_topog_refinedSampling.py --hgridfilename ~/python_tools/Analysis/mom6/TSS/Analysis/ocean_hgrid.nc  --outputfilename ocean_topog_hgrid.nc --source_file  ~/dataset/world_bathy/bathy_marmara.nc  --source_lon x --source_lat y --source_elv z

df = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/TSS/ocean_topog_hgrid_marmara_bathy.nc.nc')
lon_rho = np.copy(df['x'][1::2,1::2])
lat_rho = np.copy(df['y'][1::2,1::2])

bathy = np.copy(df.height[1::2,1::2])
bathy[np.where(bathy<0)]=0
bathy[np.where(bathy==0)]=np.nan

# northwest side of black sea coast
bathy[1461:1469,141:147] = np.nan
# minimum depth of 2 meters
bathy[np.where(bathy<2)]=2
# maximum depth of 2085 meters
bathy[np.where(bathy>2085)]=2085

# set nan to zero
bathy[np.isnan(bathy)]=0.0

nj,ni = bathy.shape
print('nj=%i, nj=%i'%(nj,ni))

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



