import numpy as np
import cartopy.crs as ccrs
import cartopy.feature
import mpl_util
from cpttoseg import cpt2seg
import matplotlib as mpllib


_palette_data = cpt2seg('/home/milicak/python_tools/Analysis/NorESM/FAMOS/Analysis/coldblue.cpt')
palette = mpllib.colors.LinearSegmentedColormap('palette', _palette_data, 11)

# df = xr.open_dataset('/archive/milicak/MITgcm_c65/Projects/Arctic_4km/Exp02_0/3DArcticOcean_monthly_THETA_2017_1-12.nc')
df = xr.open_dataset('/archive/milicak/MITgcm_c65/Projects/Arctic_4km/Exp02_0/grid.nc')
dfs = xr.open_dataset('/archive/milicak/MITgcm_c65/Projects/Arctic_4km/ncfiles/2DArcticOcean_1999_06.nc')
lon = np.copy(dfs.XC)
lat = np.copy(dfs.YC)

# contours
conts = [-5000,-4000,-3000,-2500,-2000,-1500,-1000,-500,-250,-50,0]

fig, ax1 = plt.subplots(figsize=(8,8))
m = Basemap(projection='npstere',boundinglat=48,lon_0=0,resolution='l');
m.drawcoastlines();
m.fillcontinents(color='grey');
im1 = m.contourf(lon,lat,ma.masked_where(df.Depth==0,-df.Depth),conts,cmap=mpl_util.LevelColormap(conts,cmap=palette),vmax=0, vmin=-5000,latlon=True)
cb = m.colorbar(im1,"right", size="5%", pad="2%")
cb.ax.tick_params(labelsize=14)
parallels = np.arange(40.,86,5.)
# labels = [left,right,top,bottom]
m.drawparallels(parallels,labels=[True,False,False,False],fontsize=14);
meridians = np.arange(10.,351.,20.)
m.drawmeridians(meridians,labels=[False,False,False,True],fontsize=14);

# Fram Strait
lon1 = lon[496:664,580]
lat1 = lat[496:664,580]
m.plot(lon1,lat1,'r',latlon=True)
# Davis Strait
lon1 = lon[866:1012,680]
lat1 = lat[866:1012,680]
m.plot(lon1,lat1,'r',latlon=True)
# Bering Strait
lon1 = lon[659:688,1449]
lat1 = lat[659:688,1449]
m.plot(lon1,lat1,'r',latlon=True)
# Barents Sea
aa = np.arange(384,526)
bb = np.arange(386,386+71+1)
y = xr.DataArray(bb[:-1], dims='points')
x = xr.DataArray(aa[::2], dims='points')
x2 = xr.DataArray(aa[1::2], dims='points')
lon1 = np.copy(df.XC.isel(i=x,j=y))
lat1 = np.copy(df.YC.isel(i=x,j=y))
lon1 = np.concatenate((np.array([19.3]),lon1))
lat1 = np.concatenate((np.array([69.5]),lat1))
m.plot(lon1,lat1,'r',latlon=True)

plt.savefig('paperfigs/Arctic_bathy_domain.png', bbox_inches='tight',format='png',dpi=300)



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
