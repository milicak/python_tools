import numpy as np
import glob
from mpl_toolkits.basemap import Basemap
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

ds = xr.open_dataset('moc_transports.nc')
ds2 = ds.resample(time='M').mean('time')

df = xr.open_dataset('AMOC_MOM6.nc')
amocmean = df.AMOC.mean('time')
# 26.5N is df.lat[708]
amoc26N = df.AMOC[:,:,708].max('z_l')
df_annual = amoc26N.groupby("time.year").mean("time")

fig = plt.figure(figsize=(8, 8))
ax0 = plt.subplot2grid((3, 1), (0, 0), rowspan=2)
ax1 = plt.subplot2grid((3, 1), (2, 0))

im1 = ax0.pcolormesh(amocmean.lat,-amocmean.z_l,amocmean.where(amocmean!=0),
        shading='gouraud',vmin=-5,vmax=25,rasterized=True)
# add colorbar                                                             
axpos = ax0.get_position()                                                 
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0,0.03,axpos.height]);
cbar = fig.colorbar(im1, cax=cbar_ax)                             
ax0.set_ylabel(r'Depth [m]',fontsize = 12.0)
# ax0.set_xlabel(r'Lat',fontsize = 12.0)
ax0.text(25,-6900,'Lat',fontsize=12)
# cbar.ax.tick_params(labelsize=16)                                          

ax1.plot(df.time,amoc26N,'k',label='monthly')
ax1.plot(df.time[5::12],df_annual,'r',label='annual')
ax1.plot(ds2.time[:-36],ds2.moc_mar_hc10[:-36],'b',label='RAPID')
ax1.legend(loc='lower left')  
ax1.set_ylim(5,35)
ax1.grid()
ax1.set_ylabel(r'Maximum AMOC at 26.5N [Sv]',fontsize = 12.0)
ax1.set_xlabel(r'Year',fontsize = 12.0)
ax1.text(df.time[0],32,'mean = 20.34 Sv',fontsize = 12.0)


plt.savefig('paperfigs/mean_AMOC_MOM6.png', bbox_inches='tight',format='png',dpi=300)


