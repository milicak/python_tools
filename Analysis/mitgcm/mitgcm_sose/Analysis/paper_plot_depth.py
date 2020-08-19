import cartopy.crs as ccrs
import cartopy
import matplotlib.path as mpath
import matplotlib.colors as colors
import xmitgcm

fname = 'Exp01_0s_taux.nc'
df1 = xr.open_dataset(fname)

# plt.pcolormesh(df.YG,df.Z,Trxsummean*1e-6,vmin=-10,vmax=25,cmap='jet');plt.colorbar()

df = xmitgcm.open_mdsdataset('/archive/milicak/MITgcm_c65/Projects/mitgcm_sose/Exp01_0s/',geometry='sphericalpolar',prefix='THETA')

theta = np.linspace(0, 2*np.pi, 100)
center, radius = [0.5, 0.5], 0.5
verts = np.vstack([np.sin(theta), np.cos(theta)]).T
circle = mpath.Path(verts * radius + center)



dmin = -3.0
dmax = 3.0
divnorm = colors.DivergingNorm(vmin=dmin, vcenter=0, vmax=dmax)
fig, axes = plt.subplots(figsize=(7,7))
fig.canvas.draw()
plt.tight_layout()

ax1 = plt.subplot2grid(shape=(1,1),loc=(0,0), projection=ccrs.SouthPolarStereo())

ax1.set_extent([-180, 180, -90, -24.5], ccrs.PlateCarree())
ax1.add_feature(cartopy.feature.LAND,color='grey');
ax1.coastlines(resolution='50m');
c1 = ax1.pcolormesh(df.XC,df.YC,ma.masked_where(df.Depth==0,-df.Depth),
                    vmin=-5700,vmax=0,cmap='Blues_r', #cmap='viridis',
                    transform=ccrs.PlateCarree(), rasterized=True);
ax1.set_boundary(circle, transform=ax1.transAxes)

axpos = ax1.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0+0.025,0.03,axpos.height*0.95])
cbar = fig.colorbar(c1, cax=cbar_ax)


plt.savefig('paperfigs/depth.png', bbox_inches='tight',format='png',dpi=300)

