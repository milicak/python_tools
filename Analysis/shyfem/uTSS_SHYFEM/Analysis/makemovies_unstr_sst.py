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
 
fyear = 521;
lyear = 600;

ls1 = sorted(glob.glob('/work/opa/mi19918/Projects/uTSS_SHYFEM/work/OUT/*nos*'))


# for ind in xrange(fyear,fyear+1):
plt.ioff()
for ind in range(fyear,lyear+1):
    print(ind)
    sdate = "%4.4d" % (ind)
    df = xr.open_dataset(ls1[ind])
    plt.figure(figsize=(12,8))
    m = Basemap(llcrnrlon=22.5,llcrnrlat=38.5,urcrnrlon=32.,urcrnrlat=43.5,\
                rsphere=(6378137.00,6356752.3142),\
                resolution='h',projection='merc',\
                lat_0=40.,lon_0=20.,lat_ts=20.)
    
    m.drawcoastlines(linewidth=0.2)
    m.fillcontinents(color='grey')
    m.drawparallels(np.arange(38,44,1),labels=[1,1,0,0])
    m.drawmeridians(np.arange(22,33,2),labels=[0,0,0,1])
    
    if ind == fyear:
        longitude,latitude = m(np.copy(df.longitude),np.copy(df.latitude))

    im1=plt.tripcolor(longitude,latitude,df.element_index-1,
                      df.temperature[0,:,0],cmap='nice_gfdl',vmin=10,vmax=32,shading='gouraud')
    cb = m.colorbar(im1,"right", size="5%", pad="10%",ticks=[0, 0.05, 0.1, 0.15,
                                                             0.2, 0.25, 0.3, 0.35]) # pad is the distance between colorbar and figure
    # cb.set_label('$C$',rotation=0,y=1.0,labelpad=-45)
    time=str(df.time.data)
    plt.title(time[2:12])
    printname = 'gifs/sst_'+sdate+'.png'
    plt.savefig(printname, bbox_inches='tight',format='png',dpi=300)
    plt.close()




