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

# expid = 'Exp01.2'
# expid = 'Exp_20160101'
expid = 'Exp_2016_analysis_newTSIC'

fname = root_folder+project_name+'/'+expid+'/OUT_2016/'+'uTSS_lobc_chunk_0230.nos.nc'

df = xr.open_dataset(fname)
var = df.temperature[0,:,0]

df['element_index'] -= 1

plt.figure()
m = Basemap(llcrnrlon=22.5,llcrnrlat=38.5,urcrnrlon=32.,urcrnrlat=43.5,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='h',projection='merc',\
            lat_0=40.,lon_0=20.,lat_ts=20.)

m.drawcoastlines(linewidth=0.2)
m.fillcontinents(color='grey')
m.drawparallels(np.arange(38,44,1),labels=[1,0,0,0])
m.drawmeridians(np.arange(22,33,2),labels=[0,0,0,1])

longitude,latitude = m(np.copy(df.longitude),np.copy(df.latitude))

im1=plt.tripcolor(longitude,latitude,df.element_index,
                np.ma.masked_equal(var,0),cmap='needJet2',vmin=22,vmax=28,shading='gouraud')

cb = m.colorbar(im1,"right", size="5%", pad="10%")
# cb.set_label('$C$',rotation=0,y=1.07,labelpad=-45)

# plt.savefig('paperfigs/uTSS_SSS.eps', bbox_inches='tight',format='eps', dpi=300)
plt.savefig('paperfigs/uTSS_SST_Marmara_18_08_2016.png', bbox_inches='tight',format='png', dpi=300)
plt.close()

plt.figure()
m = Basemap(llcrnrlon=26.0,llcrnrlat=40.0,urcrnrlon=30.,urcrnrlat=41.2,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='h',projection='merc',\
            lat_0=40.,lon_0=20.,lat_ts=20.)

m.drawcoastlines(linewidth=0.2)
m.fillcontinents(color='grey')
m.drawparallels(np.arange(38,44,1),labels=[1,0,0,0])
m.drawmeridians(np.arange(22,33,2),labels=[0,0,0,1])

longitude,latitude = m(np.copy(df.longitude),np.copy(df.latitude))

im1=plt.tripcolor(longitude,latitude,df.element_index,
                np.ma.masked_equal(var,0),cmap='needJet2',vmin=22,vmax=28,shading='gouraud')

cb = m.colorbar(im1,"right", size="5%", pad="10%")
# cb.set_label('$C$',rotation=0,y=1.07,labelpad=-45)
plt.savefig('paperfigs/uTSS_SST_Marmara_18_08_2016_zoom.png', bbox_inches='tight',format='png', dpi=300)
plt.close()


############# 22th August 2016 ##############
fname = root_folder+project_name+'/'+expid+'/OUT_2016/'+'uTSS_lobc_chunk_0234.nos.nc'

df = xr.open_dataset(fname)
var = df.temperature[0,:,0]

df['element_index'] -= 1

plt.figure()
m = Basemap(llcrnrlon=22.5,llcrnrlat=38.5,urcrnrlon=32.,urcrnrlat=43.5,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='h',projection='merc',\
            lat_0=40.,lon_0=20.,lat_ts=20.)

m.drawcoastlines(linewidth=0.2)
m.fillcontinents(color='grey')
m.drawparallels(np.arange(38,44,1),labels=[1,0,0,0])
m.drawmeridians(np.arange(22,33,2),labels=[0,0,0,1])

longitude,latitude = m(np.copy(df.longitude),np.copy(df.latitude))

im1=plt.tripcolor(longitude,latitude,df.element_index,
                np.ma.masked_equal(var,0),cmap='needJet2',vmin=22,vmax=28,shading='gouraud')

cb = m.colorbar(im1,"right", size="5%", pad="10%")
# cb.set_label('$C$',rotation=0,y=1.07,labelpad=-45)

# plt.savefig('paperfigs/uTSS_SSS.eps', bbox_inches='tight',format='eps', dpi=300)
plt.savefig('paperfigs/uTSS_SST_Marmara_22_08_2016.png', bbox_inches='tight',format='png', dpi=300)
plt.close()

plt.figure()
m = Basemap(llcrnrlon=26.0,llcrnrlat=40.0,urcrnrlon=30.,urcrnrlat=41.2,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='h',projection='merc',\
            lat_0=40.,lon_0=20.,lat_ts=20.)

m.drawcoastlines(linewidth=0.2)
m.fillcontinents(color='grey')
m.drawparallels(np.arange(38,44,1),labels=[1,0,0,0])
m.drawmeridians(np.arange(22,33,2),labels=[0,0,0,1])

longitude,latitude = m(np.copy(df.longitude),np.copy(df.latitude))

im1=plt.tripcolor(longitude,latitude,df.element_index,
                np.ma.masked_equal(var,0),cmap='needJet2',vmin=22,vmax=28,shading='gouraud')

cb = m.colorbar(im1,"right", size="5%", pad="10%")
# cb.set_label('$C$',rotation=0,y=1.07,labelpad=-45)
plt.savefig('paperfigs/uTSS_SST_Marmara_22_08_2016_zoom.png', bbox_inches='tight',format='png', dpi=300)


