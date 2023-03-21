import numpy as np
import cartopy.crs as ccrs
import cartopy.feature
import mpl_util
from cpttoseg import cpt2seg
import matplotlib as mpllib
import matplotlib.font_manager


_palette_data = cpt2seg('/home/milicak/python_tools/Analysis/NorESM/FAMOS/Analysis/coldblue.cpt')
palette = mpllib.colors.LinearSegmentedColormap('palette', _palette_data, 11)

df = xr.open_dataset('/archive/milicak/dataset/world_bathy/GEBCO_2014_1D.nc')
lon = np.arange(df.x_range[0],df.x_range[1],df.spacing[0])
lat = np.arange(df.y_range[0],df.y_range[1],df.spacing[1])
Depth = np.copy(df.z)
Depth.shape = np.flip(np.copy(df.dimension))
Depth[Depth>=0] = 0


lon,lat = np.meshgrid(lon,lat)

skipx = 10
skipy = 20
lon = lon[::skipx,::skipy]
lat = lat[::skipx,::skipy]
Depth = Depth[::skipx,::skipy]
Depth = np.flipud(Depth)

# contours
# conts = [-5000,-4000,-3000,-2500,-2000,-1500,-1000,-500,-250,-50,0]
conts = [-5000,-4500,-4000,-3500,-3000,-2500,-2000,-1500,-1000,-500,-400,-250,-100,-50,0]

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

# Kanada Baseni
x,y=m(-120,75)
plt.annotate('Kanada Baseni', xy=(x, y),  xycoords='data', rotation=-30,
             color='white',weight='bold')
x,y=m(-85,86)
plt.annotate('Makarov Baseni', xy=(x, y),  xycoords='data', rotation=50,
             color='white',weight='bold',fontsize=8)
x,y=m(-52,85)
# plt.annotate('Lomonosov Sırtı', xy=(x, y),  xycoords='data', rotation=43,
             # color='white',weight='bold')
plt.annotate('Lomonosov Sırtı', xy=(x, y),  xycoords='data', rotation=43,
             color='white',weight='bold',fontname='Times New Roman',fontsize=14)
# x,y=m(5,85)
# plt.annotate('Avrasya Baseni', xy=(x, y),  xycoords='data', rotation=50,
             # color='white',weight='bold')
x,y=m(-30,85)
plt.annotate('Amundsen Baseni', xy=(x, y),  xycoords='data', rotation=45,
             color='white',weight='bold',fontsize=8)
x,y=m(14,84)
plt.annotate('Nansen Baseni', xy=(x, y),  xycoords='data', rotation=55,
             color='white',weight='bold',fontsize=8)
x,y=m(21,73)
plt.annotate('Barents Denizi', xy=(x, y),  xycoords='data', rotation=-45,
             color='black',weight='bold')
x,y=m(-17,71)
plt.annotate('Nordik Denizi', xy=(x, y),  xycoords='data', rotation=0,
             color='black',weight='bold')
x,y=m(-40,70)
plt.annotate('Grönland', xy=(x, y),  xycoords='data', rotation=50,
             color='red',weight='bold',fontname='Times New Roman',fontsize=16)

x,y=m(112,78)
plt.annotate('Laptev Denizi', xy=(x, y),  xycoords='data', rotation=90,
             color='black',weight='bold',fontsize=8)
x,y=m(-164,76)
plt.annotate('Chukchi Denizi', xy=(x, y),  xycoords='data', rotation=100,
             color='black',weight='bold',fontsize=8)

# Fram Strait
x,y=m(np.array([10,-20]),np.array([79.5,79.5]))
m.plot(x,y,'r',linewidth=3)
# Bering Strait
x,y=m(np.array([-170,-167]),np.array([66,66]))
m.plot(x,y,'g',linewidth=3)
# Barents Sea
x,y=m(np.array([19.0,15.5]),np.array([69.8,77.0]))
m.plot(x,y,'k',linewidth=3)
# Davis Strait
x,y=m(np.array([-53.1,-62.0]),np.array([66.8,66.8]))
m.plot(x,y,'m',linewidth=3)


# plt.savefig('paperfigs/Arctic_real_bathy_domain.png', bbox_inches='tight',format='png',dpi=300)
plt.savefig('paperfigs/Arctic_real_bathy_domain_empty.pdf', bbox_inches='tight',format='pdf',dpi=300)
