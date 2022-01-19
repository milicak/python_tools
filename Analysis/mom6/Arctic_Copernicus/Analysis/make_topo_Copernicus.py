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
# OMtopogen/create_topog_refinedSampling.py --hgridfilename ocean_hgrid.nc  --outputfilename ocean_topog_hgrid.nc --source_file ~/dataset/world_bathy/GEBCO_2021.nc  --source_lon lon --source_lat lat --source_elv elevation

# df = xr.open_dataset('/okyanus/users/milicak/python_libs/ocean_model_topog_generator/ocean_topog_hgrid.nc')
df = xr.open_dataset('ocean_topog_hgrid.nc')
# gr = xr.open_dataset('roms_grid_Arctic_Copernicus.nc')

aa = np.copy(df.height.where(df.height<0,0))
aa = aa[1::2,::2]
aa[(aa>-5) & (aa<0)]=-5
aa[aa<-6500]=-6500
nj,ni = aa.shape

# mask = gr.mask_rho[1:-1,1:-1]
# bathy = aa*mask

mdic = {"hraw": aa}
savemat("Arctic_bathy_raw.mat", mdic)

# for Arctic setup
print('Now you have to run matlab and remove_bays_Copernicus.m')
print('load the new matfile')
tmp = loadmat('Arctic_bathy_raw_smooth.mat')
hraw = -tmp['bath4']
# great lakes closed
hraw[0,1300:1380] = 0
# bay in siberia
hraw[-1,1350:1370] = 0

# some lakes as we know
hraw[1019:1030,860:885] = 0
hraw[824:840,636:650] = 0
hraw[912:920,1769:1780] = 0
hraw[902:914,1780:1790] = 0
hraw[875:880,1842:1848] = 0
hraw[926:930,1866:1870] = 0
hraw[918:924,1937:1944] = 0
hraw[735:755,679:696] = 0
hraw[987:998,941:952] = 0
hraw[1176:1188,1485:1496] = 0
hraw[1282:1287,1380:1390] = 0
hraw[654:662,1304:1310] = 0
hraw[700:705,1450:1460] = 0
hraw[415:420,1510:1515] = 0
hraw[456:466,1483:1490] = 0
hraw[420:426,1630:1640] = 0
hraw[415:420,1640:1650] = 0
hraw[384:388,1655:1661] = 0
hraw[330:344,1656:1663] = 0
hraw[416:420,1511:1515] = 0
hraw[250:270,1200:1222] = 0
hraw[220:240,1130:1170] = 0
hraw[310:320,1040:1048] = 0
hraw[324:322,1045:1051] = 0
hraw[80:110,1610:1630] = 0
hraw[0,670:690] = 0
hraw[686:690,3:7] = 0
hraw[974:976,1595:1600] = 0
hraw[972:976,1600:1610] = 0
hraw[924:928,1660:1665] = 0
hraw[937:940,1668:1671] = 0
hraw[140:150,1440:1450] = 0
hraw[995:1010,1997:] = 0
hraw[1058:1080,1966:1976] = 0
hraw[1031:1034,1831:1834] = 0
hraw[923:936,1989:] = 0
hraw[1197:1203,1248:1251]=5
hraw[1197,1251]=5

# to plot
bb = np.where(hraw == 0, np.nan, hraw)

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

# a lake in the Canada region is closed
# also river end in the Siberia
df = xr.open_dataset('ocean_topog.nc')
# great lakes closed
hraw[0,1300:1380] = 0

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

# ncap2 -s 'defdim("ntiles",1)' ocean_topog.nc dnm.nc
cmnd1 ="ncap2 -s " + " 'defdim(" + '"ntiles"' + ",1)' " + " ocean_topog.nc dnm.nc"
os.system(cmnd1)
cmnd1 = 'mv dnm.nc ocean_topog.nc'
os.system(cmnd1)
cmnd1 = 'cp ocean_topog.nc ~/dataset/MOM6/Arctic_Copernicus/.'
os.system(cmnd1)

# this is to generatre ocean_mask.nc land_mask.nc tiles etc

cmnd1 = '/okyanus/users/milicak/models/FRE-NCtools//tools/make_quick_mosaic/make_quick_mosaic --input_mosaic ocean_mosaic.nc --mosaic_name grid_spec --ocean_topog ocean_topog.nc'


# ncap2 -s 'defdim("ntiles",1)' ocean_topog.nc dnm.nc
# mv dnm.nc ocean_topog.nc
# cp ocean_topog.nc ~/dataset/MOM6/Arctic/.

