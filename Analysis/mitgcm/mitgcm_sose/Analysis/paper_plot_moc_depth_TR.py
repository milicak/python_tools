import matplotlib.colors as colors

fname = 'Exp01_0s_MOC_depth.nc'
df1 = xr.open_dataset(fname)
fname = 'Exp02_0s_MOC_depth.nc'
df3 = xr.open_dataset(fname)

dmin = -5.0
dmax = 20.0
divnorm = colors.DivergingNorm(vmin=dmin, vcenter=0, vmax=dmax)

fig, axes = plt.subplots(figsize=(14,6))
ax1 = plt.subplot2grid(shape=(1,4),loc=(0,0), colspan=2)
ax2 = plt.subplot2grid((1,4), (0,2), colspan=2)
plt.tight_layout()

c1 = ax1.pcolormesh(df1.YG,df1.Z,df1.MOC_depth*1e-6,vmin=-10,vmax=25,cmap='needJet2');
axpos = ax1.get_position()
cbar_ax = fig.add_axes([axpos.x0-0.05,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(c1, cax=cbar_ax, ticklocation='left')
ax1.set_yticklabels([]);
ax1.set_yticks([]) ;
c2 = ax2.pcolormesh(df1.YG,df1.Z,(df3.MOC_depth-df1.MOC_depth)*1e-6,vmin=dmin,vmax=dmax,cmap='RdBu_r',norm=divnorm)
axpos = ax2.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(c2, cax=cbar_ax)
ax1.set_xlabel('Enlem', fontsize=16)
ax2.set_xlabel('Enlem', fontsize=16)
ax1.tick_params(axis='x', which='major', labelsize=16)
ax2.tick_params(axis='x', which='major', labelsize=16)
ax1.text(-77,-5400,'a)',fontsize=24)
ax2.text(-77,-5400,'b)',fontsize=24)

plt.savefig('paperfigs_tr/MOC_depth.png', bbox_inches='tight',format='png',dpi=300)

