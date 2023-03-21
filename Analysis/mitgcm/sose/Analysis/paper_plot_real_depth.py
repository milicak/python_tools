import cartopy.crs as ccrs
import cartopy
import matplotlib.path as mpath
import matplotlib.colors as colors
from matplotlib import gridspec
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


lon,lat = np.meshgrid(lon,lat)

skipx = 10
skipy = 20
lon = lon[::skipx,::skipy]
lat = lat[::skipx,::skipy]
Depth = Depth[::skipx,::skipy]
Depth = np.flipud(Depth)

# conts = [-8000, -6000, -5500, -5000,-4500,-4000,-3500,-3000,-2500,-2000,-1500,-1000,-500,-400,-250,-100,-50,0]
# fig, ax1 = plt.subplots(figsize=(8,8))
# m = Basemap(projection='spstere',boundinglat=-30,lon_0=180,resolution='h');
# m.drawcoastlines();
# m.fillcontinents(color='grey');
# # im1 = m.contourf(lon,lat,ma.masked_where(Depth==0,Depth),conts,cmap=mpl_util.LevelColormap(conts,cmap=palette),vmax=0, vmin=-5000,latlon=True)
# im1 = m.contourf(lon,lat,ma.masked_where(Depth==0,Depth),conts,cmap='Blues_r',vmax=0,
#            vmin=-6000,latlon=True)
# cb = m.colorbar(im1,"right", size="5%", pad="2%")
# cb.ax.tick_params(labelsize=14)
# parallels = np.arange(40.,86,5.)
# # labels = [left,right,top,bottom]
# m.drawparallels(parallels,labels=[True,False,False,False],fontsize=14);
# meridians = np.arange(10.,351.,20.)
# m.drawmeridians(meridians,labels=[False,False,False,True],fontsize=14);

theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)
fig = plt.figure(figsize=(8, 8))
gs = gridspec.GridSpec(1, 1)
ax0 = plt.subplot(gs[0],projection=ccrs.SouthPolarStereo())
ax0.set_extent([-180, 180, -90, -30], ccrs.PlateCarree())
ax0.add_feature(cartopy.feature.LAND,color='grey');
ax0.coastlines(resolution='50m');
c1 = ax0.pcolormesh(lon,lat,ma.masked_where(Depth==0,Depth),
                    vmin=-5700,vmax=0,cmap='Blues_r', #cmap='viridis',
                    transform=ccrs.PlateCarree(), rasterized=True);
axpos = ax0.get_position()
# cbar_ax = fig.add_axes([axpos.x0-0.05,axpos.y0,0.03,axpos.height])
cbar_ax = fig.add_axes([axpos.x1+0.02,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(c1, cax=cbar_ax)
# cbar_ax.yaxis.set_ticks_position('left')
cbar_ax.tick_params(labelsize=14)
# ax0.text(-48, -9, 'a)', transform=ccrs.Geodetic(), fontsize=14)


# plt.savefig('paperfigs/Antarctic_real_depth.png', bbox_inches='tight',format='png',dpi=300)
plt.savefig('paperfigs/Antarctic_real_depth.pdf',
            bbox_inches='tight',format='pdf',dpi=300)

