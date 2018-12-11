import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import ESMF
from mpl_toolkits.basemap import Basemap                                            
#import cartopy.crs as ccrs                                                          
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
 
plt.ion()

root_folder  = '/work/mi19918/Projects/'
project_name = 'uTSS'

# expid = 'Exp01.2'
expid = 'Exp_20160101'

aa = np.loadtxt('thalwegCoords_sorted.txt')

lon_thlweg = aa[:,0]
lat_thlweg = aa[:,1]

fyear = 0;
lyear = 30;

sdate = "%c%4.4d%c" % ('*',fyear,'*')
fname = root_folder+project_name+'/'+expid+'/OUT/'+'uTSS_lobc_chunk_'+sdate+'.nos.nc'
list = sorted(glob.glob(fname))
str1 = ''.join(list)
grd = xr.open_dataset(str1)
for year in xrange(fyear+1,lyear+1):
    sdate = "%c%4.4d%c" % ('*',year,'*')
    list.extend(sorted(glob.glob(root_folder+project_name+'/'+expid+'/OUT/'+'uTSS_lobc_chunk_'+sdate+'.nos.nc')))



data = xr.open_mfdataset(list, chunks={'time':5, 'node':2000})
data['element_index'] -= 1
grd['element_index'] -= 1

coord_sys = ESMF.CoordSys.SPH_DEG
domask = True
# create locstream
locstream = ESMF.LocStream(lon_thlweg.shape[0], name="uTSS Thalweg Section", coord_sys=coord_sys)
# appoint the section locations
locstream["ESMF:Lon"] = lon_thlweg
locstream["ESMF:Lat"] = lat_thlweg
if domask:
    locstream["ESMF:Mask"] = np.array(np.ones(lon_thlweg.shape[0]), dtype=np.int32)


secfield = np.zeros((lon_thlweg.shape[0],np.copy(data.level.shape[0])))

# triang=mtri.Triangulation(grd.longitude,grd.latitude,grd.element_index)


