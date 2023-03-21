import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import os
from datetime import date
from eofs.xarray import Eof

root_folder = '/archive/milicak/MITgcm_c65/Projects/Arctic_4km/'

expid0 = 'Exp02_0';
# Atlantic v1
expid1 = 'Exp02_1';
# Atlantic v2
expid3 = 'Exp02_3';

gridname = root_folder + expid0 + '/grid.nc'
gr = xr.open_dataset(gridname)

df1 = xr.open_dataset('/archive/milicak/MITgcm_c65/Projects/Arctic_4km/ncfiles/Exp02_1_seaice_area_diff_EOF.nc')
df3 = xr.open_dataset('/archive/milicak/MITgcm_c65/Projects/Arctic_4km/ncfiles/Exp02_3_seaice_area_diff_EOF.nc')

m = Basemap(projection='npstere',boundinglat=48,lon_0=0,resolution='l');
lon1,lat1 = m(-132,38)
lon2,lat2 = m(-45,67)
lon = np.copy(gr.XC)
lat = np.copy(gr.YC)

fig, axes = plt.subplots(figsize=(9,6))
ax1 = plt.subplot2grid(shape=(1,4),loc=(0,0),colspan=2)
ax2 = plt.subplot2grid(shape=(1,4),loc=(0,2),colspan=2)
plt.tight_layout()

m.ax=ax1
m.drawcoastlines()
m.fillcontinents(color='grey');
im1 = m.pcolormesh(lon,lat,ma.masked_where(gr.Depth==0,df1.eof1[0,:,:]),
                   cmap=plt.cm.get_cmap('RdBu_r',16),vmin=-1,vmax=1,latlon=True, rasterized=True)
parallels = np.arange(40.,86,5.)
# labels = [left,right,top,bottom]
m.drawparallels(parallels,labels=[True,False,False,False],fontsize=14);
meridians = np.arange(10.,351.,20.)
m.drawmeridians(meridians,labels=[False,False,False,True],fontsize=14);
ax1.text(lon1,lat1,'a)',fontsize=14)
ax1.text(lon2,lat2,'14%',fontsize=14,color='r')

m.ax=ax2
m.drawcoastlines()
m.fillcontinents(color='grey');
im2 = m.pcolormesh(lon,lat,ma.masked_where(gr.Depth==0,df3.eof1[0,:,:]),
                   cmap=plt.cm.get_cmap('RdBu_r',16),vmin=-1,vmax=1,latlon=True, rasterized=True)
parallels = np.arange(40.,86,5.)
# labels = [left,right,top,bottom]
m.drawparallels(parallels,labels=[False,False,False,False],fontsize=14);
meridians = np.arange(10.,351.,20.)
m.drawmeridians(meridians,labels=[False,False,False,True],fontsize=14);
ax2.text(lon1,lat1,'b)',fontsize=14)
ax2.text(lon2,lat2,'10%',fontsize=14,color='r')
axpos = ax2.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.02,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(im2, cax=cbar_ax, ticklocation='right')
cbar.ax.tick_params(labelsize=14)


plt.savefig('paperfigs/sea_ice_diff_EOF.png', bbox_inches='tight',format='png',dpi=300)
