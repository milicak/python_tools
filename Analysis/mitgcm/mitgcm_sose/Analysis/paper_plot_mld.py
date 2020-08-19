import cartopy.crs as ccrs
import xmitgcm
import cartopy
import matplotlib.path as mpath
import matplotlib.colors as colors

fname = 'Exp01_0s_taux.nc'
df1 = xr.open_dataset(fname)
fname = 'Exp02_0s_taux.nc'
df2 = xr.open_dataset(fname)
fname = 'Exp03_0s_taux.nc'
df3 = xr.open_dataset(fname)
fname = 'Exp04_0s_taux.nc'
df4 = xr.open_dataset(fname)
fname = 'Exp05_0s_taux.nc'
df5 = xr.open_dataset(fname)

# plt.pcolormesh(df.YG,df.Z,Trxsummean*1e-6,vmin=-10,vmax=25,cmap='jet');plt.colorbar()


theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)



dmin = -3.0
dmax = 3.0
divnorm = colors.DivergingNorm(vmin=dmin, vcenter=0, vmax=dmax)
fig, axes = plt.subplots(figsize=(14,12))
ax1 = plt.subplot2grid(shape=(2,6),loc=(0,1), colspan=2, projection=ccrs.SouthPolarStereo())
ax2 = plt.subplot2grid((2,6), (0,3), colspan=2, projection=ccrs.SouthPolarStereo())
ax3 = plt.subplot2grid((2,6), (1,0), colspan=2, projection=ccrs.SouthPolarStereo())
ax4 = plt.subplot2grid((2,6), (1,2), colspan=2, projection=ccrs.SouthPolarStereo())
ax5 = plt.subplot2grid((2,6), (1,4), colspan=2, projection=ccrs.SouthPolarStereo())
fig.canvas.draw()
plt.tight_layout()

ax1.set_extent([-180, 180, -90, -24.5], ccrs.PlateCarree())
ax1.add_feature(cartopy.feature.LAND,color='grey');
ax1.coastlines(resolution='50m');
c1 = ax1.pcolormesh(df1.XG,df1.YC,ma.masked_equal(df1.taux,0),
                    vmin=-.3,vmax=.3,cmap='viridis',
                    transform=ccrs.PlateCarree(), rasterized=True);
ax1.set_boundary(circle, transform=ax1.transAxes)

ax2.set_extent([-180, 180, -90, -24.5], ccrs.PlateCarree())
ax2.add_feature(cartopy.feature.LAND,color='grey');
ax2.coastlines(resolution='50m');
c2 = ax2.pcolormesh(df1.XG,df1.YC,ma.masked_equal(df2.taux,0),
                    vmin=-.3,vmax=.3,cmap='viridis',
                    transform=ccrs.PlateCarree(), rasterized=True);
ax2.set_boundary(circle, transform=ax2.transAxes)

ax3.set_extent([-180, 180, -90, -24.5], ccrs.PlateCarree())
ax3.add_feature(cartopy.feature.LAND,color='grey');
ax3.coastlines(resolution='50m');
ax3.pcolormesh(df1.XG,df1.YC,ma.masked_equal(df3.taux,0),
                    vmin=-.3,vmax=.3,cmap='viridis',
                    transform=ccrs.PlateCarree(), rasterized=True);
ax3.set_boundary(circle, transform=ax3.transAxes)

ax4.set_extent([-180, 180, -90, -24.5], ccrs.PlateCarree())
ax4.add_feature(cartopy.feature.LAND,color='grey');
ax4.coastlines(resolution='50m');
ax4.pcolormesh(df1.XG,df1.YC,ma.masked_equal(df4.taux,0),
                    vmin=-.3,vmax=.3,cmap='viridis',
                    transform=ccrs.PlateCarree(), rasterized=True);
ax4.set_boundary(circle, transform=ax4.transAxes)

ax5.set_extent([-180, 180, -90, -24.5], ccrs.PlateCarree())
ax5.add_feature(cartopy.feature.LAND,color='grey');
ax5.coastlines(resolution='50m');
ax5.pcolormesh(df1.XG,df1.YC,ma.masked_equal(df5.taux,0),
                    vmin=-.3,vmax=.3,cmap='viridis',
                    transform=ccrs.PlateCarree(), rasterized=True);
ax5.set_boundary(circle, transform=ax5.transAxes)

axpos = ax2.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0+0.025,0.03,axpos.height*0.95])
cbar = fig.colorbar(c2, cax=cbar_ax)
# axpos = ax1.get_position()
# cbar_ax = fig.add_axes([axpos.x1-0.36,axpos.y0+0.025,0.03,axpos.height*0.95])
# cbar = fig.colorbar(c1, cax=cbar_ax, ticklocation='left')

plt.savefig('paperfigs/taux.png', bbox_inches='tight',format='png',dpi=300)

