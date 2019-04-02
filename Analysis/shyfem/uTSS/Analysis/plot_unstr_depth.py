import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib import rc
from needJet2 import shfn
# import ESMF
from mpl_toolkits.basemap import Basemap                                            
#import cartopy.crs as ccrs                                                          
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
 
plt.ion()
cmap_needjet2=shfn()
#rc('text', usetex=True)

root_folder  = '/work/mi19918/Projects/'
project_name = 'uTSS'

expid = 'Exp01.2'

fname = root_folder+project_name+'/'+expid+'/'+'uTSSm0_ous.nc'

data = xr.open_dataset(fname)
data['element_index'] -= 1
#data = xr.open_dataset(fname)['salinity']

triang = mtri.Triangulation(data.longitude,data.latitude,data.element_index)


plt.figure()
m = Basemap(llcrnrlon=22.5,llcrnrlat=38.5,urcrnrlon=32.,urcrnrlat=43.5,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='h',projection='merc',\
            lat_0=40.,lon_0=20.,lat_ts=20.)

m.drawcoastlines(linewidth=0.2)
m.fillcontinents(color='grey')
m.drawparallels(np.arange(38,44,1),labels=[1,1,0,0])
m.drawmeridians(np.arange(22,33,2),labels=[0,0,0,1])

longitude,latitude = m(np.copy(data.longitude),np.copy(data.latitude))

im1=plt.tripcolor(longitude,latitude,data.element_index,
                -data.total_depth,cmap=cmap_needjet2,vmin=-2100,vmax=-2,shading='gouraud')
#plt.tripcolor(longitude,latitude,data.element_index,ksinode/fzero,cmap=plt.cm.get_cmap('RdBu_r'),vmin=-1,vmax=1,shading='gouraud')


cb = m.colorbar(im1,"right", size="5%", pad="10%",
                extend='min')
                                                 
cb.set_label('$[m]$',rotation=0,y=1.07,labelpad=-45)
plt.savefig('paperfigs/uTSS_bathy.png', bbox_inches='tight',format='png',dpi=300)

plt.close()

plt.figure()
m = Basemap(llcrnrlon=28.8,llcrnrlat=40.9,urcrnrlon=29.3,urcrnrlat=41.3,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='h',projection='merc',\
            lat_0=40.,lon_0=20.,lat_ts=20.)

m.drawcoastlines(linewidth=0.2)
m.fillcontinents(color='grey')
m.drawparallels(np.arange(40.8,41.4,.2),labels=[1,1,0,0])
m.drawmeridians(np.arange(28.8,29.4,.25),labels=[0,0,0,1])


longitude,latitude = m(np.copy(data.longitude),np.copy(data.latitude))
im1=plt.tripcolor(longitude,latitude,data.element_index,
                -data.total_depth,cmap=cmap_needjet2,vmin=-100,vmax=-2,shading='gouraud')


cb = m.colorbar(im1,"right", size="5%", pad="10%",
                extend='min')
                                                 
cb.set_label('$[m]$',rotation=0,y=1.07,labelpad=-45)
plt.savefig('paperfigs/uTSS_bosphorus_bathy.png', bbox_inches='tight',format='png',dpi=300)
