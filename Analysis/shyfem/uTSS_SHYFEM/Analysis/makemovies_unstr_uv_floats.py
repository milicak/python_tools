import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib import rc
# import ESMF
from mpl_toolkits.basemap import Basemap                                            
from datetime import date, timedelta
#import cartopy.crs as ccrs                                                          
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib.colors import LinearSegmentedColormap

colors = np.loadtxt('/users_home/opa/mi19918/.python/nice_gfdl.rgb')                     
cm = LinearSegmentedColormap.from_list('my_list', colors, N=225)   
plt.register_cmap('nice_gfdl',cm)                                               
 
fyear = 0;
lyear = 143;

df = xr.open_dataset('/work/opa/mi19918/Projects/uTSS_SHYFEM/work/out_0729/dnm2.nc')
ds = xr.open_dataset('/work/opa/mi19918/Projects/uTSS_SHYFEM/work/out_0729/opendrift.nc')
ds2 = xr.open_dataset('/work/opa/mi19918/Projects/uTSS_SHYFEM/work/out_0729/opendrift_hormix.nc')
# compute surface speed
sp = np.copy(df.u_velocity[:,:,0])**2 + np.copy(df.v_velocity[:,:,0])**2
sp = np.sqrt(sp)

# for ind in xrange(fyear,fyear+1):
plt.ioff()
m = Basemap(llcrnrlon=26,llcrnrlat=40.0,urcrnrlon=30.,urcrnrlat=41.5,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='h',projection='merc',\
            lat_0=40.,lon_0=20.,lat_ts=20.)
m.drawcoastlines(linewidth=0.2)
m.fillcontinents(color='grey')
m.drawparallels(np.arange(38,44,1),labels=[1,1,0,0])
m.drawmeridians(np.arange(22,33,2),labels=[0,0,0,1])
for ind in range(fyear,lyear+1):
    print(ind)
    sdate = "%4.4d" % (ind)
    plt.figure(figsize=(9,6))
    if ind == fyear:
        longitude,latitude = m(np.copy(df.longitude),np.copy(df.latitude))

    im1=plt.tripcolor(longitude,latitude,df.element_index-1,
                     sp[round(ind/2),:],cmap='nice_gfdl',vmin=0,vmax=1,shading='gouraud')
    cb = m.colorbar(im1,"right", size="5%", pad="10%")
    m.scatter(ds.lon[:,ind],ds.lat[:,ind],c='k',latlon=True,s=3,label='No diffusion')
    m.scatter(ds2.lon[:,ind],ds2.lat[:,ind],c='r',latlon=True,s=3,label='Diffusion')
    plt.legend(loc='upper left')
    # cb.set_label('$C$',rotation=0,y=1.0,labelpad=-45)
    time=str(ds.time[ind].data)
    plt.title(time[2:16])
    printname = 'gifs/uvfloats_'+sdate+'.png'
    plt.savefig(printname, bbox_inches='tight',format='png',dpi=300)
    plt.close()




