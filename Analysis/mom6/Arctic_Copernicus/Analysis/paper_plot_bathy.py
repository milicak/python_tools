import numpy as np
import cartopy.crs as ccrs
import cartopy.feature
import mpl_util
from cpttoseg import cpt2seg
import matplotlib as mpllib
import cmocean


_palette_data = cpt2seg('/home/milicak/python_tools/Analysis/NorESM/FAMOS/Analysis/coldblue.cpt')
palette = mpllib.colors.LinearSegmentedColormap('palette', _palette_data, 11)

dfs = xr.open_dataset('/archive/milicak/MOM6-examples/Projects/Arctic_Copernicus/ocean_geometry.nc')
lon = np.copy(dfs.geolon)
lat = np.copy(dfs.geolat)

# contours
# conts = [-6500, -5000,-4000,-3000,-2500,-2000,-1500,-1000,-500,-250,-50,0]
conts = [-6500,-5000,-4000,-3000,-2500,-2000,-1500,-1000,-500,-400,-350,-300,-250,-200,-150,-100,-50,0]

fig, ax1 = plt.subplots(figsize=(8,8))
m = Basemap(projection='npstere',boundinglat=37,lon_0=0,resolution='l');
m.drawcoastlines();
m.fillcontinents(color='grey');
im1 = m.contourf(lon,lat,ma.masked_where(dfs.D==0,-dfs.D),conts,cmap=mpl_util.LevelColormap(conts,cmap=cmocean.cm.deep_r),vmax=0, vmin=-5000,latlon=True)
# im1 = m.contourf(lon,lat,ma.masked_where(dfs.D==0,-dfs.D),conts,cmap=mpl_util.LevelColormap(conts,cmap=palette),vmax=0, vmin=-5000,latlon=True)
cb = m.colorbar(im1,"right", size="5%", pad="2%")
cb.ax.tick_params(labelsize=14)
parallels = np.arange(40.,86,5.)
# labels = [left,right,top,bottom]
m.drawparallels(parallels,labels=[True,False,False,False],fontsize=14);
meridians = np.arange(10.,351.,20.)
m.drawmeridians(meridians,labels=[False,False,False,True],fontsize=14);

# Fram Strait
# Fram Strait indices are j=724:844, i = 674
lon1 = lon[724:844,1424]
lat1 = lat[724:844,1424]
m.plot(lon1,lat1,'k',latlon=True,linewidth=2)
# Bering Strait
# Bering Strait indices are j=689:713, i = 674
lon1 = lon[689:713,674]
lat1 = lat[689:713,674]
m.plot(lon1,lat1,'r',latlon=True,linewidth=2)
# Davis Strait
aa = np.arange(1450,1501)
bb = np.arange(350,401)
y = xr.DataArray(bb, dims='points')
x = xr.DataArray(aa, dims='points')
lon1 = np.copy(dfs.geolon.isel(lonh=x,lath=y))
lat1 = np.copy(dfs.geolat.isel(lonh=x,lath=y))
m.plot(lon1,lat1,'g',latlon=True,linewidth=2)
# Barents Sea
# aa = np.arange(384,526)
# bb = np.arange(386,386+71+1)
# y = xr.DataArray(bb[:-1], dims='points')
# x = xr.DataArray(aa[::2], dims='points')
# x2 = xr.DataArray(aa[1::2], dims='points')
# lon1 = np.copy(df.XC.isel(i=x,j=y))
# lat1 = np.copy(df.YC.isel(i=x,j=y))
# lon1 = np.concatenate((np.array([19.3]),lon1))
# lat1 = np.concatenate((np.array([69.5]),lat1))
# m.plot(lon1,lat1,'r',latlon=True)
#
plt.savefig('paperfigs/Arctic_bathy_domain.png', bbox_inches='tight',format='png',dpi=300)

