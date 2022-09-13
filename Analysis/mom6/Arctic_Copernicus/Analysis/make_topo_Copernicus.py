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
hraw[0,1300:1380] = -0
# bay in siberia
hraw[-1,1350:1370] = -0

# some lakes as we know
hraw[1019:1030,860:885] = -0
hraw[824:840,636:650] = -0
hraw[912:920,1769:1780] = -0
hraw[902:914,1780:1790] = -0
hraw[875:880,1842:1848] = -0
hraw[926:930,1866:1870] = -0
hraw[918:924,1937:1944] = -0
hraw[735:755,679:696] = -0
hraw[987:998,941:952] = -0
hraw[1176:1188,1485:1496] = -0
hraw[1282:1287,1380:1390] = -0
hraw[654:662,1304:1310] = -0
hraw[700:705,1450:1460] = -0
hraw[415:420,1510:1515] = -0
hraw[456:466,1483:1490] = -0
hraw[420:426,1630:1640] = -0
hraw[415:420,1640:1650] = -0
hraw[384:388,1655:1661] = -0
hraw[330:344,1656:1663] = -0
hraw[416:420,1511:1515] = -0
hraw[250:270,1200:1222] = -0
hraw[220:240,1130:1170] = -0
hraw[310:320,1040:1048] = -0
hraw[324:322,1045:1051] = -0
hraw[80:110,1610:1630] = -0
hraw[0,670:690] = -0
hraw[686:690,3:7] = -0
hraw[974:976,1595:1600] = -0
hraw[972:976,1600:1610] = -0
hraw[924:928,1660:1665] = -0
hraw[937:940,1668:1671] = -0
hraw[140:150,1440:1450] = -0
hraw[995:1010,1997:] = -0
hraw[1058:1080,1966:1976] = -0
hraw[1031:1034,1831:1834] = -0
hraw[923:936,1989:] = -0
hraw[1197:1203,1248:1251]=5
hraw[1197,1251]=5
hraw[735:760,1996:]=-0
hraw[705:709,1988:1990]=-0
hraw[0:20,1740:1799]=-0
hraw[0:2,1625:1635]=-0
hraw[3:6,687:690]=-0
hraw[85:90,667:670]=-0
hraw[134:136,670:673]=-0
hraw[144:146,673:676]=-0
hraw[199:203,647:651]=-0
hraw[206:209,700:704]=-0
hraw[624:634,707:714]=-0
hraw[637:639,660:663]=-0
hraw[670:674,687:690]=-0
hraw[646:649,593:595]=-0
hraw[567:569,453:455]=-0
hraw[550:552,475:477]=-0
hraw[607:609,407:409]=-0
hraw[706:708,1988:1990]=-0
hraw[1074:1076,895:897]=-0
hraw[1016:1018,479:481]=-0
hraw[344:350,1067:1073]=-0
hraw[388:392,1236:1240]=-0
hraw[524:526,1242:1244]=-0
hraw[587:589,1263:1265]=-0
hraw[462:467,887:903]=-0
hraw[325:332,1046:1052]=-0
hraw[362:364,1563:1565]=-0
hraw[314:322,1134:1142]=-0
hraw[358:360,1256]=hraw[358:360,1257]
hraw[358:360,1255]=hraw[358:360,1257]
hraw[358:360,1260]=hraw[358:360,1259]
hraw[1138:1140,1221:1223]=-0
hraw[1073:1075,1631:1633]=-0
hraw[1008:1010,1896:1898]=-0
hraw[472:485,1622:1628]=-0
hraw[476:478,1454:1456]=-0
hraw[508:510,1401:1403]=-0
hraw[1016:1018,1183:1185]=-0
hraw[868:870,1817:1819]=-0
hraw[692:694,1481:1483]=-0
hraw[473:477,1628:1630]=-0
hraw[175:177,669:671]=-0
hraw[860:862,1847:1849]=-0
hraw[628:630,816:818]=-0
hraw[1287:1294,1377:1379]=-0


# weird east boundary island at the boundary
hraw[686:690,3]=0.5*(hraw[691,3]+hraw[685,3])
hraw[686:690,4]=0.5*(hraw[691,4]+hraw[685,4])
hraw[686:690,5]=0.5*(hraw[691,5]+hraw[685,5])
hraw[686:690,6]=0.5*(hraw[691,6]+hraw[685,6])

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
hraw = np.copy(df.depth)
# great lakes closed
hraw[0,1300:1380] = 0

# 4 points lake closed
hraw[161:163,675:677]=-0.0
nj,ni = hraw.shape

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

