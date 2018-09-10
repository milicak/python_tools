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

fname = root_folder+project_name+'/'+expid+'/'+'uTSSm0_ous.nc'

data = xr.open_dataset(fname)
data['element_index'] -= 1
#data = xr.open_dataset(fname)['salinity']

triang = mtri.Triangulation(data.longitude,data.latitude,data.element_index)

# velocity time,node,level indexing
varu = data.u_velocity[-1,:,0]
varv = data.v_velocity[-1,:,0]

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

plt.tripcolor(longitude,latitude,data.element_index,
                varu,cmap='jet',vmin=-1,vmax=1,shading='gouraud')

#plt.savefig('paperfigs/uTSS_uvel.png', bbox_inches='tight',format='png',dpi=300)
