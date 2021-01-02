import numpy as np
import cartopy.crs as ccrs
import cartopy.feature
import mpl_util
from cpttoseg import cpt2seg
import matplotlib as mpllib


dfs = xr.open_dataset('/archive/milicak/MITgcm_c65/Projects/Arctic_4km/ncfiles/2DArcticOcean_1999_06.nc')
lon = np.copy(dfs.XC)
lat = np.copy(dfs.YC)
df = xr.open_dataset('/archive/milicak/MITgcm_c65/Projects/Arctic_4km/ncfiles/Exp02_0_SST.nc')
dnm = df.sel(year=slice(2007,2016))
sst = dnm.mean('year')

fig = plt.figure(figsize=(6,6))
fig.canvas.draw()
plt.tight_layout()
m = Basemap(projection='npstere',boundinglat=48,lon_0=0,resolution='l');
# x,y = m(-45, 34.5)
m.drawcoastlines();
m.fillcontinents(color='grey');
im1 = m.pcolormesh(lon,lat,ma.masked_where(df.Depth==0,sst.sst),cmap='needJet2',vmin=-2,
           vmax=12,latlon=True, rasterized=True)
cb = m.colorbar(im1,"right", size="5%", pad="4%")
cb.ax.tick_params(labelsize=14)
parallels = np.arange(40.,86,5.)
# labels = [left,right,top,bottom]
m.drawparallels(parallels,labels=[True,False,False,False],fontsize=14);
meridians = np.arange(10.,351.,20.)
m.drawmeridians(meridians,labels=[False,False,False,True],fontsize=14);

plt.savefig('paperfigs/ctrl_sst_mean.png', bbox_inches='tight',format='png',dpi=300)

# m = Basemap(projection='npstere',boundinglat=48,lon_0=0,resolution='l',
#             ax=axes[0]);




# m = Basemap(width=6000000,height=6000000,
#          resolution='l',projection='stere',
# 	 lat_ts=45,lat_0=90,lon_0=0.,round='true')
# m.drawparallels(parallels)
# m.drawmeridians(meridians)
# m.drawcoastlines()
# m.fillcontinents(color='1')
# m.drawrivers()
#
#
#
# fig, ax1 = plt.subplots(figsize=(8,8))
# ax1 = plt.subplot(1, 1, 1, projection=ccrs.NorthPolarStereo())
# ax1.set_extent([-180, 180, 50, 90], ccrs.PlateCarree())
# ax1.add_feature(cartopy.feature.LAND,color='grey');
# ax1.coastlines(resolution='50m');
# c1 = ax1.pcolormesh(df.XC,df.YC,ma.masked_where(df.Depth==0,-df.Depth),
#                     # cmap=mpl_util.LevelColormap(conts,cmap=palette),vmax=0, vmin=-5000,
#                     cmap='needJet2',vmax=0, vmin=-5000,
#                     transform=ccrs.PlateCarree(), rasterized=True);
#
# # fig.colorbar(c1, orientation="horizontal", pad=0.1)
# axpos = ax1.get_position()
# cbar_ax = fig.add_axes([axpos.x1+0.1,axpos.y0,0.03,axpos.height*0.95])
# cbar = fig.colorbar(c1, cax=axpos)
#
#
# c1 = ax1.contourf(df.XC,df.YC,ma.masked_where(df.Depth==0,-df.Depth),
#                     conts,cmap=mpl_util.LevelColormap(conts,cmap=palette),vmax=0, vmin=-5000,
#                     transform=ccrs.PlateCarree());
#
#
