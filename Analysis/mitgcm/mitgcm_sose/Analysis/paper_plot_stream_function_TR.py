import cartopy.crs as ccrs
import cartopy
import matplotlib.path as mpath
import matplotlib.colors as colors

fname = 'Exp01_0s_stream_function.nc'
df1 = xr.open_dataset(fname)
fname = 'Exp03_0s_stream_function.nc'
df3 = xr.open_dataset(fname)

theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)


dmin = -5.0
dmax = 20.0
divnorm = colors.DivergingNorm(vmin=dmin, vcenter=0, vmax=dmax)

fig, axes = plt.subplots(figsize=(14,12))
ax1 = plt.subplot2grid(shape=(1,4),loc=(0,0), colspan=2, projection=ccrs.SouthPolarStereo())
ax2 = plt.subplot2grid((1,4), (0,2), colspan=2, projection=ccrs.SouthPolarStereo())
fig.canvas.draw()
plt.tight_layout()

ax1.set_extent([-180, 180, -90, -24.5], ccrs.PlateCarree())
ax1.add_feature(cartopy.feature.LAND,color='grey');
ax1.coastlines(resolution='50m');
c1 = ax1.pcolormesh(df1.XG,df1.YC,ma.masked_equal(df1.stream_function*1e-6,0),
                    vmin=-100,vmax=250,cmap='needJet2',
                    transform=ccrs.PlateCarree(), rasterized=True);
ax1.set_boundary(circle, transform=ax1.transAxes)

ax2.set_extent([-180, 180, -90, -24.5], ccrs.PlateCarree())
ax2.add_feature(cartopy.feature.LAND,color='grey');
ax2.coastlines(resolution='50m');
c2 = ax2.pcolormesh(df1.XG,df1.YC,ma.masked_equal(df3.stream_function*1e-6,0)
                    -ma.masked_equal(df1.stream_function*1e-6,0),
                    vmin=-50,vmax=50,cmap='RdBu_r', #norm=divnorm,
                    transform=ccrs.PlateCarree(), rasterized=True);
ax2.set_boundary(circle, transform=ax2.transAxes)


axpos = ax2.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0+0.025,0.03,axpos.height*0.95])
cbar = fig.colorbar(c2, cax=cbar_ax)
axpos = ax1.get_position()
cbar_ax = fig.add_axes([axpos.x0-0.05,axpos.y0+0.025,0.03,axpos.height*0.95])
cbar = fig.colorbar(c1, cax=cbar_ax, ticklocation='left')
ax1.text(-40, -15, 'a)', transform=ccrs.Geodetic(), fontsize=30);
ax2.text(-40, -15, 'b)', transform=ccrs.Geodetic(), fontsize=30);

plt.savefig('paperfigs_tr/stream_function.png', bbox_inches='tight',format='png',dpi=300)

