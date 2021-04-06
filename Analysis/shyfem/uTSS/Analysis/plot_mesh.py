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

m = Basemap(llcrnrlon=22.5,llcrnrlat=38.5,urcrnrlon=32.,urcrnrlat=43.5,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='h',projection='merc',\
            lat_0=40.,lon_0=20.,lat_ts=20.)

m.drawcoastlines(linewidth=0.2);
m.fillcontinents(color='grey');
m.drawparallels(np.arange(38,44,1),labels=[1,1,0,0]);
m.drawmeridians(np.arange(22,33,2),labels=[0,0,0,1]);

longitude,latitude = m(np.copy(df.longitude[0,:]),np.copy(df.latitude[0,:]))

im1 = plt.tripcolor(longitude,latitude,df.element_index[-1,:,:]-1,np.sqrt(sp)
             ,vmin=0,vmax=0.6,shading='gouraud',cmap=cmocean.cm.speed);

cb = m.colorbar(im1,"right", size="5%", pad="10%", extend='max')
cb.set_label('$[m/s]$',rotation=0,y=1.07,labelpad=-45)
plt.savefig('paperfigs/uTSS_bathy_mesh.png', bbox_inches='tight',format='png',dpi=300)


