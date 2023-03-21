import numpy as np
import cmocean

df1 = xr.open_dataset('Canada_TS_mom6.nc')
df2 = xr.open_dataset('Canada_TS_WOA13.nc')

ds = df1.temp_Canada[264-10*12:].mean('time')

f, (a0, a1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [3, 1]})

a1.plot(ds, -ds.depth,'k',label='MOM6')
a1.plot(df2.temp_Canada,-df2.depth,'r',label='WOA18')
a1.set_ylim(-4100,0)
a1.set_xlabel(r'$^\circ$C',fontsize = 14.0)
a1.set_yticklabels([])
# a1.grid()
a1.legend()

im2 = a0.pcolormesh(df1.time,-df1.depth,np.transpose(np.copy(df1.temp_Canada)),
                  cmap=cmocean.cm.thermal)
a0.set_ylim(-4100,0)
a0.set_xlabel(r'year',fontsize = 14.0)
f.tight_layout()

axpos = a0.get_position()
cbar_ax = f.add_axes([axpos.x0-0.14,axpos.y0,0.03,axpos.height])
cbar = f.colorbar(im2, cax=cbar_ax, ticklocation='left')
cbar.ax.tick_params(labelsize=14)
cbar.set_label('$^\circ$C',rotation=0,y=1.02,labelpad=-40,fontsize=14)
# a0.text(lon1,lat1,'a)',fontsize=14);

f.savefig('paperfigs/Canada_temp.png', bbox_inches='tight',format='png',dpi=300)

