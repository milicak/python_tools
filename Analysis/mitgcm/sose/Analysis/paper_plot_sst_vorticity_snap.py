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

def coriolis(lat):
    """Compute the Coriolis parameter for the given latitude:
    ``f = 2*omega*sin(lat)``, where omega is the angular velocity
    of the Earth.

    Parameters
    ----------
    lat : array
      Latitude [degrees].
    """
    omega   = 7.2921159e-05  # angular velocity of the Earth [rad/s]
    return 2*omega*np.sin(lat/360.*2*np.pi)


# plt figure
df1 = xr.open_dataset('/archive2/milicak/mitgcm/sose/Exp01_0/SST_2011_10_25.nc')
df2 = xr.open_dataset('/archive2/milicak/mitgcm/sose/Exp01_0/vort2dsurf_2011_10_25.nc')
gr = xr.open_dataset('/archive2/milicak/mitgcm/sose/Exp01_0/grid.nc')
fcor = coriolis(np.copy(gr.YC))
aa = np.matlib.repmat(fcor, 4320, 1)
fcor = np.transpose(aa)

fig, axes = plt.subplots(figsize=(9,6))
ax1 = plt.subplot2grid(shape=(1,4),loc=(0,0),
                       colspan=2,projection=ccrs.SouthPolarStereo())
ax2 = plt.subplot2grid(shape=(1,4),loc=(0,2),
                       colspan=2,projection=ccrs.SouthPolarStereo())

ax1.set_extent([-180, 180, -90, -29.5], ccrs.PlateCarree())
ax1.add_feature(cartopy.feature.LAND,color='grey');
ax1.coastlines(resolution='110m');
ax2.set_extent([-180, 180, -90, -29.5], ccrs.PlateCarree())
ax2.add_feature(cartopy.feature.LAND,color='grey');
ax2.coastlines(resolution='110m');
# fig.canvas.draw()
# wait first to plot the figure then the use command below
# plt.tight_layout()

c1 = ax1.pcolormesh(gr.XC,gr.YC,ma.masked_where(gr.Depth==0,df1.THETA[0,:,:]),
                    vmin=-2,vmax=25,cmap='nice_gfdl',
                    transform=ccrs.PlateCarree(), rasterized=True);
axpos = ax1.get_position()
cbar_ax = fig.add_axes([axpos.x0-0.05,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(c1, cax=cbar_ax, ticklocation='left')
cbar.ax.tick_params(labelsize=14)
cbar.set_label('$^\circ$C',rotation=0,y=1.01,labelpad=-31,fontsize=14)
ax1.text(-50, -18, 'a)', transform=ccrs.Geodetic(), fontsize=14);

c2 = ax2.pcolormesh(gr.XC,gr.YC,ma.masked_where(gr.Depth==0,
                                                df2.momVort3[0,:,:])/fcor,
                    vmin=-0.25,vmax=0.25,cmap='BrBG',
                    transform=ccrs.PlateCarree(), rasterized=True);
axpos = ax2.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.03,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(c2, cax=cbar_ax, ticklocation='right')
cbar.ax.tick_params(labelsize=14)
cbar.set_label(r'f/$\zeta$',rotation=0,y=1.07,labelpad=-48,fontsize=14)
ax2.text(-50, -18, 'b)', transform=ccrs.Geodetic(), fontsize=14);

plt.savefig('paperfigs/snap_SST_vort.png', bbox_inches='tight',format='png',dpi=300)

