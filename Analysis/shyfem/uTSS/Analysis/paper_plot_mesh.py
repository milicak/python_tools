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
import cmocean
import matplotlib.colors as colors
#import cartopy.crs as ccrs
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

dmin = -2.0
dmax = 10.0
# divnorm = colors.TwoSlopeNorm(vmin=dmin, vcenter=0, vmax=dmax)
divnorm = colors.DivergingNorm(vmin=dmin, vcenter=0, vmax=dmax)

aa = np.loadtxt('newthalweg.txt')

fname = '/archive/milicak/shyfem/uTSS/Exp_2016_analysis_newTSIC/OUT/uTSS_lobc_chunk_1400.nos.nc'
df = xr.open_dataset(fname)

fig = plt.figure(figsize=(10, 6))
ax1 = fig.add_axes([0.05, 0.05, 0.95, 0.95])
m = Basemap(llcrnrlon=22.5,llcrnrlat=38.5,urcrnrlon=32.,urcrnrlat=43.5,
            rsphere=(6378137.00,6356752.3142),
            resolution='h',projection='merc',
            lat_0=40.,lon_0=20.,lat_ts=20.)

m.drawcoastlines(linewidth=0.2);
m.fillcontinents(color='grey');
m.drawparallels(np.arange(38,44,1),labels=[1,1,0,0]);
m.drawmeridians(np.arange(22,33,2),labels=[0,0,0,1]);
longitude,latitude = m(np.copy(df.longitude),np.copy(df.latitude))
im1 = plt.tripcolor(longitude,latitude,df.element_index-1,df.total_depth
             ,vmin=0,vmax=1200,edgecolor='k',linewidth=0.05,cmap=cmocean.cm.deep);
lon1,lat1 = m(aa[:,0],aa[:,1])
plt.plot(lon1,lat1,'r');

cb = m.colorbar(im1,"right", size="5%", pad="10%", extend='max')
cb.set_label('[m]',rotation=0,y=1.03,labelpad=-45)

# Bozburun and Saroz
x, y = m(26.5, 40.5)
x2, y2 = (-30, 30.5)
plt.annotate('Saroz Bay', xy=(x, y),  xycoords='data',
                xytext=(x2, y2), textcoords='offset points',
                color='k',
                arrowprops=dict(arrowstyle="fancy", color='g')
                )

x, y = m(28.8, 40.53)
x2, y2 = (90,30.5)
plt.annotate('Bozburun Peninsula', xy=(x, y),  xycoords='data',
                xytext=(x2, y2), textcoords='offset points',
                color='k',
                arrowprops=dict(arrowstyle="fancy", color='g')
                )

x, y = m(28.53, 40.52)
x2, y2 = (-100,-50.95)
plt.annotate('Imrali Island', xy=(x, y),  xycoords='data',
                xytext=(x2, y2), textcoords='offset points',
                color='k',
                arrowprops=dict(arrowstyle="fancy", color='g')
                )

ax2 = fig.add_axes([0.05, 0.6, 0.37, 0.37])
m = Basemap(llcrnrlon=25.8,llcrnrlat=39.8,urcrnrlon=27.,urcrnrlat=40.6,
            rsphere=(6378137.00,6356752.3142),
            resolution='h',projection='merc',
            lat_0=40.,lon_0=20.,lat_ts=20.)
m.drawcoastlines(linewidth=0.2);
m.drawparallels(np.arange(38,44,0.5),labels=[0,1,0,0]);
m.drawmeridians(np.arange(22,33,0.5),labels=[0,0,0,1]);
longitude,latitude = m(np.copy(df.longitude),np.copy(df.latitude))
im1 = plt.tripcolor(longitude,latitude,df.element_index-1,df.total_depth
             ,vmin=0,vmax=100,shading='gouraud',cmap='RdYlBu_r');
cb = m.colorbar(im1,"right", size="5%", pad="20%", extend='max')
cb.set_label('[m]',rotation=0,y=1.07,labelpad=-29)

# Southern Section at Dardanelles
Nx = 50
lat_SC_Dar = np.linspace(40,40.045,Nx)
lon_SC_Dar = 26.188*np.ones(lat_SC_Dar.shape)
x, y = m(lon_SC_Dar,lat_SC_Dar)
plt.plot(x,y,'k',linewidth=4)
# Northern Section at Dardanelles
lon_NC_Dar = np.array([26.6979,26.7375])
lat_NC_Dar = np.array([40.4208,40.3937])
x, y = m(lon_NC_Dar,lat_NC_Dar)
plt.plot(x,y,'k',linewidth=4)


ax3 = fig.add_axes([0.55, 0.05, 0.35, 0.35])
m = Basemap(llcrnrlon=28.8,llcrnrlat=40.9,urcrnrlon=29.3,urcrnrlat=41.3,
            rsphere=(6378137.00,6356752.3142),
            resolution='h',projection='merc',
            lat_0=40.,lon_0=20.,lat_ts=20.)
m.drawcoastlines(linewidth=0.2);
m.drawparallels(np.arange(38,44,0.2),labels=[1,0,0,0]);
m.drawmeridians(np.arange(22,33,0.2),labels=[0,0,1,0]);
longitude,latitude = m(np.copy(df.longitude),np.copy(df.latitude))
im1 = plt.tripcolor(longitude,latitude,df.element_index-1,df.total_depth
             ,vmin=0,vmax=80,shading='gouraud',cmap='RdYlBu_r');
cb = m.colorbar(im1,"left", size="5%", pad="30%", extend='max')
cb.set_label('[m]',rotation=0,y=0.98,labelpad=9)
cb.ax.yaxis.set_ticks_position("left")

# Southern Section at Bosphorus
lon_SC_Bos = np.linspace(28.984,29.013,Nx)
lat_SC_Bos = 41.0123*np.ones(lon_SC_Bos.shape)
x, y = m(lon_SC_Bos,lat_SC_Bos)
plt.plot(x,y,'k',linewidth=4)

lon_NC_Bos = np.linspace(29.102,29.142,Nx)
lat_NC_Bos = 41.212*np.ones(lon_NC_Bos.shape)
x, y = m(lon_NC_Bos,lat_NC_Bos)
plt.plot(x,y,'k',linewidth=4)


plt.savefig('paperfigs/uTSS_bathy_mesh.png', bbox_inches='tight',format='png',dpi=300)

