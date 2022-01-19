# This file is designed to be cut and pasted into an ipython --pylab
# session. Otherwise, you'll need to "import np as np" then
# convert "array" to "np.array".
import os
import numpy as np
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



# first run the following

# cd ~/python_libs/ocean_model_topog_generator
# OMtopogen/create_topog_refinedSampling.py --hgridfilename
# ~/dataset/MOM6/Arctic/ocean_hgrid.nc  --outputfilename ocean_topog_hgrid.nc
# --source_file ~/dataset/world_bathy/GEBCO_2021.nc  --source_lon lon
# --source_lat lat --source_elv elevation

# df = xr.open_dataset('/okyanus/users/milicak/python_libs/ocean_model_topog_generator/ocean_topog_hgrid.nc')
df = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic/ocean_topog_hgrid.nc')

aa = np.copy(df.height.where(df.height<0,0))
aa = aa[1::2,::2]
aa[(aa>-5) & (aa<0)]=-5
aa[aa<-6000]=-6000
nj,ni = aa.shape

mdic = {"hraw": aa}
savemat("Arctic_bathy_raw.mat", mdic)

# for Arctic setup
print('Now you have to run matlab and ')
print('load the new matfile')
tmp = loadmat('Arctic_bathy_raw_smooth.mat')
hraw = -tmp['bath4']

# After that run the code there are still blow ups you can use the code below
# grep -Eo 'i=.{0,5}' Arctic_remove2 > Arctic_i_ind
# grep -Eo 'j=.{0,5}' Arctic_remove2 > Arctic_j_ind
# cut -c3- Arctic_i_ind > Arctic_i_ind2


iind = np.loadtxt("Arctic_i_ind01",dtype='i')
jind = np.loadtxt("Arctic_j_ind01",dtype='i')
for itr in range(0,len(iind)):
    hraw[jind[itr]-1,iind[itr]-1] = -0.0

iind = np.loadtxt("Arctic_i_ind02",dtype='i')
jind = np.loadtxt("Arctic_j_ind02",dtype='i')
for itr in range(0,len(iind)):
    hraw[jind[itr]-1,iind[itr]-1] = -0.0

iind = np.loadtxt("Arctic_i_ind03",dtype='i')
jind = np.loadtxt("Arctic_j_ind03",dtype='i')
for itr in range(0,len(iind)):
    hraw[jind[itr]-1,iind[itr]-1] = -0.0

iind = np.loadtxt("Arctic_i_ind04",dtype='i')
jind = np.loadtxt("Arctic_j_ind04",dtype='i')
for itr in range(0,len(iind)):
    hraw[jind[itr]-1,iind[itr]-1] = -0.0

iind = np.loadtxt("Arctic_i_ind05",dtype='i')
jind = np.loadtxt("Arctic_j_ind05",dtype='i')
for itr in range(0,len(iind)):
    hraw[jind[itr]-1,iind[itr]-1] = -0.0

iind = np.loadtxt("Arctic_i_ind06",dtype='i')
jind = np.loadtxt("Arctic_j_ind06",dtype='i')
for itr in range(0,len(iind)):
    hraw[jind[itr]-1,iind[itr]-1] = -0.0

iind = np.loadtxt("Arctic_i_ind07",dtype='i')
jind = np.loadtxt("Arctic_j_ind07",dtype='i')
for itr in range(0,len(iind)):
    hraw[jind[itr]-1,iind[itr]-1] = -0.0

iind = np.loadtxt("Arctic_i_ind08",dtype='i')
jind = np.loadtxt("Arctic_j_ind08",dtype='i')
for itr in range(0,len(iind)):
    hraw[jind[itr]-1,iind[itr]-1] = -0.0

iind = np.loadtxt("Arctic_i_ind09",dtype='i')
jind = np.loadtxt("Arctic_j_ind09",dtype='i')
for itr in range(0,len(iind)):
    hraw[jind[itr]-1,iind[itr]-1] = -0.0

iind = np.loadtxt("Arctic_i_ind10",dtype='i')
jind = np.loadtxt("Arctic_j_ind10",dtype='i')
for itr in range(0,len(iind)):
    hraw[jind[itr]-1,iind[itr]-1] = -0.0

iind = np.loadtxt("Arctic_i_ind11",dtype='i')
jind = np.loadtxt("Arctic_j_ind11",dtype='i')
for itr in range(0,len(iind)):
    hraw[jind[itr]-1,iind[itr]-1] = -0.0

iind = np.loadtxt("Arctic_i_ind12",dtype='i')
jind = np.loadtxt("Arctic_j_ind12",dtype='i')
for itr in range(0,len(iind)):
    hraw[jind[itr]-1,iind[itr]-1] = -0.0

iind = np.loadtxt("Arctic_i_ind13",dtype='i')
jind = np.loadtxt("Arctic_j_ind13",dtype='i')
for itr in range(0,len(iind)):
    hraw[jind[itr]-1,iind[itr]-1] = -0.0

iind = np.loadtxt("Arctic_i_ind14",dtype='i')
jind = np.loadtxt("Arctic_j_ind14",dtype='i')
for itr in range(0,len(iind)):
    hraw[jind[itr]-1,iind[itr]-1] = -0.0

iind = np.loadtxt("Arctic_i_ind15",dtype='i')
jind = np.loadtxt("Arctic_j_ind15",dtype='i')
for itr in range(0,len(iind)):
    hraw[jind[itr]-1,iind[itr]-1] = -0.0

iind = np.loadtxt("Arctic_i_ind16",dtype='i')
jind = np.loadtxt("Arctic_j_ind16",dtype='i')
for itr in range(0,len(iind)):
    hraw[jind[itr]-1,iind[itr]-1] = -0.0

iind = np.loadtxt("Arctic_i_ind17",dtype='i')
jind = np.loadtxt("Arctic_j_ind17",dtype='i')
for itr in range(0,len(iind)):
    hraw[jind[itr]-1,iind[itr]-1] = -0.0

iind = np.loadtxt("Arctic_i_ind18",dtype='i')
jind = np.loadtxt("Arctic_j_ind18",dtype='i')
for itr in range(0,len(iind)):
    hraw[jind[itr]-1,iind[itr]-1] = -0.0

iind = np.loadtxt("Arctic_i_ind19",dtype='i')
jind = np.loadtxt("Arctic_j_ind19",dtype='i')
for itr in range(0,len(iind)):
    hraw[jind[itr]-1,iind[itr]-1] = -0.0

iind = np.loadtxt("Arctic_i_ind20",dtype='i')
jind = np.loadtxt("Arctic_j_ind20",dtype='i')
for itr in range(0,len(iind)):
    hraw[jind[itr]-1,iind[itr]-1] = -0.0

iind = np.loadtxt("Arctic_i_ind21",dtype='i')
jind = np.loadtxt("Arctic_j_ind21",dtype='i')
for itr in range(0,len(iind)):
    hraw[jind[itr]-1,iind[itr]-1] = -0.0

iind = np.loadtxt("Arctic_i_ind22",dtype='i')
jind = np.loadtxt("Arctic_j_ind22",dtype='i')
for itr in range(0,len(iind)):
    hraw[jind[itr]-1,iind[itr]-1] = -0.0

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


cmnd1 ='ncap2 -s ' + 'defdim("ntiles",1)' + 'ocean_topog.nc dnm.nc'

cmnd1 = 'mv dnm.nc ocean_topog.nc'

cmnd1 = 'mv ocean_topog.nc ~/dataset/MOM6/Arctic/.'

# ncap2 -s 'defdim("ntiles",1)' ocean_topog.nc dnm.nc
# mv dnm.nc ocean_topog.nc
# mv ocean_topog.nc ~/dataset/MOM6/Arctic/.

