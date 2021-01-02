import cartopy.crs as ccrs
import cartopy
import matplotlib.path as mpath
import matplotlib.colors as colors
from matplotlib import gridspec


df = xr.open_dataset('/archive2/milicak/mitgcm/sose/Exp01_0/grid.nc')

theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)

dmin = -3.0
dmax = 3.0
divnorm = colors.DivergingNorm(vmin=dmin, vcenter=0, vmax=dmax)


# fig, axes = plt.subplots(figsize=(7,7))
# fig.canvas.draw()
# plt.tight_layout()
# ax0 = plt.subplot2grid(shape=(1,1),loc=(0,0), projection=ccrs.SouthPolarStereo())
# ax0.set_extent([-180, 180, -90, -24.5], ccrs.PlateCarree())

fig = plt.figure(figsize=(9, 6))
gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
ax0 = plt.subplot(gs[0],projection=ccrs.SouthPolarStereo())
ax0.set_extent([-180, 180, -90, -19.5], ccrs.PlateCarree())
ax0.add_feature(cartopy.feature.LAND,color='grey');
ax0.coastlines(resolution='50m');
c1 = ax0.pcolormesh(df.XC,df.YC,ma.masked_where(df.Depth==0,-df.Depth),
                    vmin=-5700,vmax=0,cmap='Blues_r', #cmap='viridis',
                    transform=ccrs.PlateCarree(), rasterized=True);
# ax0.set_boundary(circle, transform=ax0.transAxes)
# cb.set_label('m',rotation=0,y=1.07,labelpad=-36)
axpos = ax0.get_position()
cbar_ax = fig.add_axes([axpos.x1-0.56,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(c1, cax=cbar_ax)
cbar_ax.yaxis.set_ticks_position('left')
cbar_ax.tick_params(labelsize=14)

# second plot
ax1 = plt.subplot(gs[1])
df1 = xr.open_dataset('Exp01_0_zonal_mean_stress.nc')
df2 = xr.open_dataset('Exp02_0_zonal_mean_stress.nc')
df3 = xr.open_dataset('Exp03_0_zonal_mean_stress.nc')

ax1.plot(df1.oceTAUX,df.YC,'b',label='CTRL')
ax1.plot(df3.oceTAUX,df.YC,'r',label='THERMAL')
ax1.tick_params(labelsize=14)
ax1.yaxis.tick_right()
ax1.legend(fontsize=12);
ax1.legend(frameon=False)
ax1.set_xlabel('Nm$^{-2}$',fontsize=14)
ax1.set_ylabel('Lat',fontsize=14)
ax1.yaxis.set_label_position("right")
ax1.grid()
ax0.text(-48, -9, 'a)', transform=ccrs.Geodetic(), fontsize=14)
ax1.text(-0.05,-80,'b)',fontsize=14);


plt.savefig('paperfigs/SO_depth_and_taux.png', bbox_inches='tight',format='png',dpi=300)

