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

# OMtopogen/create_topog_refinedSampling.py --hgridfilename ~/dataset/MOM6/TSS/ocean_hgrid.nc  --outputfilename ocean_topog_hgrid.nc --source_file  ~/dataset/world_bathy/gebco_East_MED_BS.nc  --source_lon lon --source_lat lat --source_elv elevation


df1 = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/TSS/ocean_topog_hgrid_marmara_bathy.nc')
df2 = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/TSS/ocean_topog_hgrid_gebco_bathy.nc')
ds = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/TSS/ocean_topog_marmara.nc')

bathy = -np.copy(df2.height[1::2,1::2])
bathy[np.where(bathy<0)]=0
bathy[np.where(bathy==0)]=np.nan

# minimum depth of 2 meters
bathy[np.where(bathy<2)]=2
# maximum depth of 2085 meters
bathy[np.where(bathy>2085)]=2085

# set nan to zero
bathy[np.isnan(bathy)]=0.0

# read the old bathy
bathy2 = np.copy(ds.depth)
# update the marmara basin
for j in range(1370,2026):
    for i in range(1375,1525):
        if bathy2[j,i] > 150:
            bathy2[j,i] = bathy[j,i]

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
hx[:,:] = bathy2
rg.close()



