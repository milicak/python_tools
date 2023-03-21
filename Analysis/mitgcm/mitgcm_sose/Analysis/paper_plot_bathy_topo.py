import cartopy.crs as ccrs
import cartopy
import matplotlib.path as mpath
import matplotlib.colors as colors


colors1 = np.loadtxt('/home/milicak/.python/OceanLakeLandSnow.rgb')
colors1 = colors1[2:,:]
colors1 /= 256
colors2 = np.loadtxt('/home/milicak/.python/GMT_ocean.rgb')

colors_t = np.vstack((colors2[:42,:], colors1))
cm1 = LinearSegmentedColormap.from_list('my_list', colors_t, N=294)
colors_t = np.vstack((colors2, colors1))
cm2 = LinearSegmentedColormap.from_list('my_list', colors_t, N=254+80)

level1=np.linspace(-10300,-62,42)
level2=np.linspace(0,7800,254)
ll=np.concatenate((level1,level2))
divnorm2 = BoundaryNorm(ll, ncolors=294, clip=True)

level1=np.linspace(-6000,-62,80)
level2=np.linspace(0,7800,254)
ll=np.concatenate((level1,level2))
divnorm3 = BoundaryNorm(ll, ncolors=254+80, clip=True)

df = xr.open_dataset('/archive/milicak/dataset/world_bathy/GEBCO_2014_1D.nc')

lon = np.linspace(df.x_range[0],df.x_range[1],np.copy(df.dimension[0]))
lat = np.linspace(df.y_range[0],df.y_range[1],np.copy(df.dimension[1]))

Depth = np.copy(df.z)
aa = np.flipud(np.reshape(Depth,(21600, 43200)))


theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)

dmin = -6000.0
dmax = 7833.0
divnorm = colors.DivergingNorm(vmin=dmin, vcenter=6000, vmax=dmax)


lon_mit = np.linspace(-180,180,360)
lat_mit = -26.5*np.ones(lon_mit.shape[0])

fig, axes = plt.subplots(figsize=(7,7))
fig.canvas.draw()
plt.tight_layout()

ax1 = plt.subplot2grid(shape=(1,1),loc=(0,0), projection=ccrs.SouthPolarStereo())
ax1.set_extent([-180, 180, -90, -14.5], ccrs.PlateCarree())
ax1.coastlines(resolution='50m');
c1 = ax1.pcolormesh(lon[::10],lat[::10],aa[::10,::10],
                    norm=divnorm2,shading='gouraud',cmap=cm1,
                    transform=ccrs.PlateCarree(), rasterized=True);
ax1.plot(lon_mit,lat_mit,'r',transform=ccrs.PlateCarree())

axpos = ax1.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0+0.025,0.03,axpos.height*0.95])
cbar = fig.colorbar(c1, cax=cbar_ax)


plt.savefig('paperfigs/bathy_topo.png', bbox_inches='tight',format='png',dpi=300)

