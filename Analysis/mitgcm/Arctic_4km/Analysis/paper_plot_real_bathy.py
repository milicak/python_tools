import numpy as np
import cartopy.crs as ccrs
import cartopy.feature
import mpl_util
from cpttoseg import cpt2seg
import matplotlib as mpllib


_palette_data = cpt2seg('/home/milicak/python_tools/Analysis/NorESM/FAMOS/Analysis/coldblue.cpt')
palette = mpllib.colors.LinearSegmentedColormap('palette', _palette_data, 11)

df = xr.open_dataset('/archive/milicak/dataset/world_bathy/GEBCO_2014_1D.nc')
lon = np.arange(df.x_range[0],df.x_range[1],df.spacing[0])
lat = np.arange(df.y_range[0],df.y_range[1],df.spacing[1])
Depth = np.copy(df.z)
Depth.shape = np.flip(np.copy(df.dimension))
Depth[Depth>=0] = 0

# contours
conts = [-5000,-4000,-3000,-2500,-2000,-1500,-1000,-500,-250,-50,0]

fig, ax1 = plt.subplots(figsize=(8,8))
m = Basemap(projection='npstere',boundinglat=60,lon_0=0,resolution='l');
m.drawcoastlines();
m.fillcontinents(color='grey');
im1 = m.contourf(lon,lat,ma.masked_where(Depth==0,Depth),conts,cmap=mpl_util.LevelColormap(conts,cmap=palette),vmax=0, vmin=-5000,latlon=True)
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
m.plot(lon1,lat1,'g',latlon=True)
# Bering Strait
lon1 = lon[659:688,1449]
lat1 = lat[659:688,1449]
m.plot(lon1,lat1,'k',latlon=True)
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
m.plot(lon1,lat1,'m',latlon=True)

plt.savefig('paperfigs/Arctic_real_bathy_domain.png', bbox_inches='tight',format='png',dpi=300)

