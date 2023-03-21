import os
import numpy as np
import cartopy.crs as ccrs
import cartopy
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib
import seaborn as sns

def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet.
        N: number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    if type(cmap) == str:
        cmap = get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in np.arange(N+1) ]
    # Return colormap object.
    return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)

# plt figure
df1=xr.open_dataset('/archive2/milicak/mitgcm/sose/Exp01_0/Exp01_0_stream_function.nc')
df2=xr.open_dataset('/archive2/milicak/mitgcm/sose/Exp02_0/Exp02_0_stream_function.nc')
df3=xr.open_dataset('/archive2/milicak/mitgcm/sose/Exp03_0/Exp03_0_stream_function.nc')
gr = xr.open_dataset('/archive2/milicak/mitgcm/sose/Exp01_0/grid.nc')

# cmap = plt.get_cmap('seismic',20)
# cmap = plt.get_cmap('RdBu_r',21)
cmap_tmp = sns.color_palette("vlag", as_cmap=True)
cmap = cmap_discretize(cmap_tmp, 15)

fig, axes = plt.subplots(figsize=(8,9))
ax1 = plt.subplot2grid(shape=(2,4),loc=(0,1),
                       colspan=2,projection=ccrs.SouthPolarStereo())
ax2 = plt.subplot2grid(shape=(2,4),loc=(1,0),
                       colspan=2,projection=ccrs.SouthPolarStereo())
ax3 = plt.subplot2grid(shape=(2,4),loc=(1,2),
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

# wait first to plot the figure then the use command below
plt.tight_layout()

c1 = ax1.pcolormesh(gr.XC,gr.YC,ma.masked_where(gr.Depth==0,df1.stream_function)*1e-6,
                    vmin=-60,vmax=200,cmap='needJet2',
                    transform=ccrs.PlateCarree(), rasterized=True);
axpos = ax1.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.03,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(c1, cax=cbar_ax, ticklocation='right')
cbar.ax.tick_params(labelsize=14)
cbar.set_label('Sv',rotation=0,y=1.07,labelpad=-44,fontsize=14)
ax1.text(-50, -18, 'a)', transform=ccrs.Geodetic(), fontsize=14);

c2 = ax2.pcolormesh(gr.XC,gr.YC,ma.masked_where(gr.Depth==0,
                                            (df2.stream_function-df1.stream_function))*1e-6,
                    vmin=-20,vmax=20,cmap=cmap,
                    transform=ccrs.PlateCarree(), rasterized=True);
ax2.text(-50, -18, 'b)', transform=ccrs.Geodetic(), fontsize=14);

c3 = ax3.pcolormesh(gr.XC,gr.YC,ma.masked_where(gr.Depth==0,
                                            (df3.stream_function-df1.stream_function))*1e-6,
                    vmin=-20,vmax=20,cmap=cmap,
                    transform=ccrs.PlateCarree(), rasterized=True);
axpos = ax3.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.03,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(c3, cax=cbar_ax, ticklocation='right')
cbar.ax.tick_params(labelsize=14)
cbar.set_label('Sv',rotation=0,y=1.07,labelpad=-44,fontsize=14)
ax3.text(-50, -18, 'c)', transform=ccrs.Geodetic(), fontsize=14);

plt.savefig('paperfigs/stream_function.png', bbox_inches='tight',format='png',dpi=300)

