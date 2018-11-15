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
expid = 'Exp_20160101'

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


triang=mtri.Triangulation(grd.longitude,grd.latitude,grd.element_index)


plt.figure()
m = Basemap(llcrnrlon=22.5,llcrnrlat=38.5,urcrnrlon=32.,urcrnrlat=43.5,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='h',projection='merc',\
            lat_0=40.,lon_0=20.,lat_ts=20.)

m.drawcoastlines(linewidth=0.2)
m.fillcontinents(color='grey')
m.drawparallels(np.arange(38,44,1),labels=[1,1,0,0])
m.drawmeridians(np.arange(22,33,2),labels=[0,0,0,1])

longitude,latitude = m(np.copy(grd.longitude),np.copy(grd.latitude))

# salinity time,node,level indexing
var = data.salinity[-1,:,0]

#plt.tricontourf(longitude,latitude,data.element_index,
#                np.ma.masked_equal(var,0),cmap='jet')
                #np.ma.masked_equal(var,0),cmap='jet',levels=range(10,40,1))
im1=plt.tripcolor(longitude,latitude,grd.element_index,
                np.ma.masked_equal(var,0),cmap='needJet2',shading='gouraud')

cb = m.colorbar(im1,"right", size="5%", pad="10%")
cb.set_label('$psu$',rotation=0,y=1.07,labelpad=-45)

# plt.savefig('paperfigs/uTSS_SSS.eps', bbox_inches='tight',format='eps', dpi=300)
plt.savefig('paperfigs/uTSS_SSS.png', bbox_inches='tight',format='png', dpi=300)
