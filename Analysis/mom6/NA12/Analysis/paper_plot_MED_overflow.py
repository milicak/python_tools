import numpy as np
import glob
from mpl_toolkits.basemap import Basemap
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

df = xr.open_dataset('MED_section_38_5_MOM6.nc')
dfs = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/WOA13/woa13_decav_s00_04v2.nc',decode_times=False)

fig = plt.figure(figsize=(8, 8))
ax0 = plt.subplot2grid((4, 1), (0, 0), rowspan=2)
ax1 = plt.subplot2grid((4, 1), (2, 0), rowspan=2)

im1 = ax0.pcolormesh(df.xh,-df.z_l,df.salt,
        vmin=34,vmax=36.5,rasterized=True)
# add colorbar                                                             
axpos = ax0.get_position()                                                 
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0,0.03,axpos.height]);
cbar = fig.colorbar(im1, cax=cbar_ax)                             
ax0.set_xlim(-25,-7.5)
ax0.set_ylabel(r'Depth [m]',fontsize = 12.0)
ax0.set_ylim(-6000,0)
# ax0.set_xlabel(r'Lat',fontsize = 12.0)
# ax0.text(25,-6900,'Lat',fontsize=12)
# cbar.ax.tick_params(labelsize=16)                                          

im1 = ax1.pcolormesh(dfs.lon,-dfs.depth,dfs.s_an[0,:,514,:],
        vmin=34,vmax=36.5,rasterized=True)
# add colorbar                                                             
axpos = ax1.get_position()                                                 
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0,0.03,axpos.height]);
cbar = fig.colorbar(im1, cax=cbar_ax)                             
ax1.set_xlim(-25,-7.5)
ax1.set_ylim(-6000,0)
ax1.set_ylabel(r'Depth [m]',fontsize = 12.0)
ax1.set_xlabel(r'Lat',fontsize = 12.0)

plt.savefig('paperfigs/MED_overflow_section_MOM6_obs.png', bbox_inches='tight',format='png',dpi=300)


