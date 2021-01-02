import matplotlib.colors as colors

fname = 'Exp01_0_MOC_depth.nc'
df1 = xr.open_dataset(fname)
fname = 'Exp02_0_MOC_depth.nc'
df2 = xr.open_dataset(fname)
fname = 'Exp03_0_MOC_depth.nc'
df3 = xr.open_dataset(fname)

# plt.pcolormesh(df.YG,df.Z,Trxsummean*1e-6,vmin=-10,vmax=25,cmap='jet');plt.colorbar()

dmin = -2.0
dmax = 18.0
depths = np.array([-6000., -5000., -4000., -3000., -2000., -1000.,     0.])
# divnorm = colors.TwoSlopeNorm(vmin=dmin, vcenter=0, vmax=dmax)
divnorm = colors.DivergingNorm(vmin=dmin, vcenter=0, vmax=dmax)

fig, axes = plt.subplots(figsize=(14,6))
ax1 = plt.subplot2grid(shape=(2,4),loc=(0,1), colspan=2)
ax2 = plt.subplot2grid((2,4), (1,0), colspan=2)
ax3 = plt.subplot2grid((2,4), (1,2), colspan=2)
plt.tight_layout()
c1 = ax1.pcolormesh(df1.YG,df1.Z,df1.MOC_depth*1e-6,vmin=-25,vmax=25,cmap='needJet2',shading='auto');
ax1.set_xlim(-78,-20);
# ax1.set_xlim(-78,-20);
axpos = ax1.get_position()
# cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0,0.03,axpos.height])
cbar_ax = fig.add_axes([axpos.x1+0.03,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(c1, cax=cbar_ax, ticklocation='right')
ax1.tick_params(labelsize=14)
ax1.set_ylabel('Depth [m]',fontsize=14)
ax1.text(-77,-5700,'a)',fontsize=14)
cbar.ax.tick_params(labelsize=14)
ax1.set_yticks([-6000,-5000,-4000,-3000,-2000,-1000,0]);
cbar.set_label('Sv',rotation=0,y=1.09,labelpad=-50,fontsize=14)


# ax1.set_yticklabels([]);
# ax1.set_yticks([]) ;
c2 = ax2.pcolormesh(df1.YG,df1.Z,(df2.MOC_depth-df1.MOC_depth)*1e-6,cmap='RdBu_r',norm=divnorm, shading='auto')
ax2.set_xlim(-78,-20);
ax2.tick_params(labelsize=14)
ax2.text(-77,-5700,'b)',fontsize=14)
ax2.set_yticks([-6000,-5000,-4000,-3000,-2000,-1000,0]);
ax2.set_xlabel('Lat',fontsize=14)
ax2.set_ylabel('Depth [m]',fontsize=14)

c3 = ax3.pcolormesh(df1.YG,df1.Z,(df3.MOC_depth-df1.MOC_depth)*1e-6,cmap='RdBu_r',norm=divnorm, shading='auto')
ax3.set_xlim(-78,-20);
axpos = ax3.get_position()
ax3.tick_params(labelsize=14)
ax3.text(-77,-5700,'c)',fontsize=14)
ax3.set_yticks([-6000,-5000,-4000,-3000,-2000,-1000,0]);
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(c3, cax=cbar_ax)
ax3.set_yticklabels([]);
ax3.set_xlabel('Lat',fontsize=14)
cbar.ax.tick_params(labelsize=14)
cbar.set_label('Sv',rotation=0,y=1.09,labelpad=-55,fontsize=14)

plt.savefig('paperfigs/MOC_depth.png', bbox_inches='tight',format='png',dpi=300)

