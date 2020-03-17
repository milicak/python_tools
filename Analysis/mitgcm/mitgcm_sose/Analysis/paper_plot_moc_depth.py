import matplotlib.colors as colors

fname = 'Exp01_0s_MOC_depth.nc'
df1 = xr.open_dataset(fname)
fname = 'Exp02_0s_MOC_depth.nc'
df2 = xr.open_dataset(fname)
fname = 'Exp03_0s_MOC_depth.nc'
df3 = xr.open_dataset(fname)
fname = 'Exp04_0s_MOC_depth.nc'
df4 = xr.open_dataset(fname)
fname = 'Exp05_0s_MOC_depth.nc'
df5 = xr.open_dataset(fname)

# plt.pcolormesh(df.YG,df.Z,Trxsummean*1e-6,vmin=-10,vmax=25,cmap='jet');plt.colorbar()

dmin = -5.0
dmax = 20.0
divnorm = colors.DivergingNorm(vmin=dmin, vcenter=0, vmax=dmax)
fig, axes = plt.subplots(figsize=(14,6))
ax1 = plt.subplot2grid(shape=(2,6),loc=(0,1), colspan=2)
ax2 = plt.subplot2grid((2,6), (0,3), colspan=2)
ax3 = plt.subplot2grid((2,6), (1,0), colspan=2)
ax4 = plt.subplot2grid((2,6), (1,2), colspan=2)
ax5 = plt.subplot2grid((2,6), (1,4), colspan=2)
plt.tight_layout()
c1 = ax1.pcolormesh(df1.YG,df1.Z,df1.MOC_depth*1e-6,vmin=-10,vmax=25,cmap='needJet2');
axpos = ax1.get_position()
# cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0,0.03,axpos.height])
cbar_ax = fig.add_axes([axpos.x1-0.33,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(c1, cax=cbar_ax, ticklocation='left')
ax1.set_yticklabels([]);
ax1.set_yticks([]) ;
c2 = ax2.pcolormesh(df1.YG,df1.Z,(df2.MOC_depth-df1.MOC_depth)*1e-6,vmin=dmin,vmax=dmax,cmap='RdBu_r',norm=divnorm)
axpos = ax2.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(c2, cax=cbar_ax)
ax3.pcolormesh(df1.YG,df1.Z,(df3.MOC_depth-df1.MOC_depth)*1e-6,vmin=dmin,vmax=dmax,cmap='RdBu_r',norm=divnorm)
ax4.pcolormesh(df1.YG,df1.Z,(df4.MOC_depth-df1.MOC_depth)*1e-6,vmin=dmin,vmax=dmax,cmap='RdBu_r',norm=divnorm)
ax5.pcolormesh(df1.YG,df1.Z,(df5.MOC_depth-df1.MOC_depth)*1e-6,vmin=dmin,vmax=dmax,cmap='RdBu_r',norm=divnorm)

plt.savefig('paperfigs/MOC_depth.png', bbox_inches='tight',format='png',dpi=300)

