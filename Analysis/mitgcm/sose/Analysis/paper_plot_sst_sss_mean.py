import os
import cartopy.crs as ccrs
import cartopy
import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import cmocean


# plt figure
df1 = xr.open_dataset('/archive2/milicak/mitgcm/sose/Exp01_0/Exp01_0_SST_SSS_mean.nc')
df2 = xr.open_dataset('/archive2/milicak/mitgcm/sose/Exp02_0/Exp02_0_SST_SSS_mean.nc')
df3 = xr.open_dataset('/archive2/milicak/mitgcm/sose/Exp03_0/Exp03_0_SST_SSS_mean.nc')
gr = xr.open_dataset('/archive2/milicak/mitgcm/sose/Exp01_0/grid.nc')

cmap = plt.get_cmap('RdBu_r',10)

fig, axes = plt.subplots(figsize=(8,15))
ax1 = plt.subplot2grid(shape=(3,4),loc=(0,0),
                       colspan=2,projection=ccrs.SouthPolarStereo())
ax2 = plt.subplot2grid(shape=(3,4),loc=(0,2),
                       colspan=2,projection=ccrs.SouthPolarStereo())
ax3 = plt.subplot2grid(shape=(3,4),loc=(1,0),
                       colspan=2,projection=ccrs.SouthPolarStereo())
ax4 = plt.subplot2grid(shape=(3,4),loc=(1,2),
                       colspan=2,projection=ccrs.SouthPolarStereo())
ax5 = plt.subplot2grid(shape=(3,4),loc=(2,0),
                       colspan=2,projection=ccrs.SouthPolarStereo())
ax6 = plt.subplot2grid(shape=(3,4),loc=(2,2),
                       colspan=2,projection=ccrs.SouthPolarStereo())

ax1.set_extent([-180, 180, -90, -29.5], ccrs.PlateCarree())
ax1.add_feature(cartopy.feature.LAND,color='grey');
ax1.coastlines(resolution='110m');
ax2.set_extent([-180, 180, -90, -29.5], ccrs.PlateCarree())
ax2.add_feature(cartopy.feature.LAND,color='grey');
ax2.coastlines(resolution='110m');
ax3.set_extent([-180, 180, -90, -29.5], ccrs.PlateCarree())
ax3.add_feature(cartopy.feature.LAND,color='grey');
ax3.coastlines(resolution='110m');
ax4.set_extent([-180, 180, -90, -29.5], ccrs.PlateCarree())
ax4.add_feature(cartopy.feature.LAND,color='grey');
ax4.coastlines(resolution='110m');
ax5.set_extent([-180, 180, -90, -29.5], ccrs.PlateCarree())
ax5.add_feature(cartopy.feature.LAND,color='grey');
ax5.coastlines(resolution='110m');
ax6.set_extent([-180, 180, -90, -29.5], ccrs.PlateCarree())
ax6.add_feature(cartopy.feature.LAND,color='grey');
ax6.coastlines(resolution='110m');
# fig.canvas.draw()
# wait first to plot the figure then the use command below
plt.tight_layout()

c1 = ax1.pcolormesh(gr.XC,gr.YC,ma.masked_where(gr.Depth==0,df1.THETA),
                    vmin=-2,vmax=25,cmap='nice_gfdl',
                    transform=ccrs.PlateCarree(), rasterized=True);
axpos = ax1.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.03,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(c1, cax=cbar_ax, ticklocation='right')
cbar.ax.tick_params(labelsize=14)
cbar.set_label('$^\circ$C',rotation=0,y=1.07,labelpad=-38,fontsize=14)
ax1.text(-50, -18, 'a)', transform=ccrs.Geodetic(), fontsize=14);

c2 = ax2.pcolormesh(gr.XC,gr.YC,ma.masked_where(gr.Depth==0,df1.SALT),
                    vmin=32.5,vmax=36.5,cmap='needJet2',
                    transform=ccrs.PlateCarree(), rasterized=True);
axpos = ax2.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.03,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(c2, cax=cbar_ax, ticklocation='right')
cbar.ax.tick_params(labelsize=14)
cbar.set_label('psu',rotation=0,y=1.07,labelpad=-33,fontsize=14)
ax2.text(-50, -18, 'b)', transform=ccrs.Geodetic(), fontsize=14);

c3 = ax3.pcolormesh(gr.XC,gr.YC,ma.masked_where(gr.Depth==0,df2.THETA-df1.THETA),
                    vmin=-1.0,vmax=1.0,cmap=cmap,
                    transform=ccrs.PlateCarree(), rasterized=True);
axpos = ax3.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.03,axpos.y0-0.14,0.03,axpos.height])
cbar = fig.colorbar(c3, cax=cbar_ax, ticklocation='right')
cbar.ax.tick_params(labelsize=14)
cbar.set_label('$^\circ$C',rotation=0,y=1.07,labelpad=-54,fontsize=14)
ax3.text(-50, -18, 'c)', transform=ccrs.Geodetic(), fontsize=14);

c4 = ax4.pcolormesh(gr.XC,gr.YC,ma.masked_where(gr.Depth==0,df2.SALT-df1.SALT),
                    vmin=-0.5,vmax=0.5,cmap=cmap,
                    transform=ccrs.PlateCarree(), rasterized=True);
axpos = ax4.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.03,axpos.y0-0.14,0.03,axpos.height])
cbar = fig.colorbar(c4, cax=cbar_ax, ticklocation='right')
cbar.ax.tick_params(labelsize=14)
cbar.set_label('psu',rotation=0,y=1.07,labelpad=-50,fontsize=14)
ax4.text(-50, -18, 'd)', transform=ccrs.Geodetic(), fontsize=14);

c5 = ax5.pcolormesh(gr.XC,gr.YC,ma.masked_where(gr.Depth==0,df3.THETA-df1.THETA),
                    vmin=-1.0,vmax=1.0,cmap=cmap,
                    transform=ccrs.PlateCarree(), rasterized=True);
ax5.text(-50, -18, 'e)', transform=ccrs.Geodetic(), fontsize=14);

c6 = ax6.pcolormesh(gr.XC,gr.YC,ma.masked_where(gr.Depth==0,df3.SALT-df1.SALT),
                    vmin=-0.5,vmax=0.5,cmap=cmap,
                    transform=ccrs.PlateCarree(), rasterized=True);
ax6.text(-50, -18, 'f)', transform=ccrs.Geodetic(), fontsize=14);


plt.savefig('paperfigs/mean_SST_SSS.png', bbox_inches='tight',format='png',dpi=300)

