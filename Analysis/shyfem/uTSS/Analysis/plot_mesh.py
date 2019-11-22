import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
# import ESMF
from mpl_toolkits.basemap import Basemap                                            
#import cartopy.crs as ccrs                                                          
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
 
plt.ion()

root_folder  = '/work/mi19918/Projects/'
project_name = 'uTSS'

expid = 'Exp01.2'
expid = 'Exp_2016_analysis'

fname = root_folder+project_name+'/'+expid+'/'+'uTSSm0_nos.nc'

data = xr.open_dataset(fname)
#data = xr.open_dataset(fname)['salinity']

triang=mtri.Triangulation(data.longitude,data.latitude,data.element_index-1)


plt.figure()
#m = Basemap(llcrnrlon=23.,llcrnrlat=38.,urcrnrlon=34.,urcrnrlat=45.,\
#            rsphere=(6378137.00,6356752.3142),\
#            resolution='h',projection='merc',\
#            lat_0=40.,lon_0=20.,lat_ts=20.)

#m.drawcoastlines(linewidth=0.2)
#m.fillcontinents(color='grey')


# salinity time,node,level indexing
plt.triplot(triang, 'k-',linewidth=0.15)
plt.savefig('paperfigs/uTSS_mesh.eps', bbox_inches='tight',format='eps', dpi=300)
plt.savefig('paperfigs/uTSS_mesh.png', bbox_inches='tight',format='png', dpi=300)
