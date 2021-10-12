import numpy as np
from mpl_toolkits.basemap import Basemap
import xarray as xr
import cmocean



gr = xr.open_dataset('/work/opa/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/OUT/uTSS_lobc_chunk_1460.ous.nc')

df1 = xr.open_dataset('/work/opa/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/daily_clim/MKE_daily.nc')
df2 = xr.open_dataset('/work/opa/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/daily_clim/EKE_daily.nc')

MKE = df1.MKE[:,275:].mean('time')
EKE = df2.EKE[275:,:].mean('time')

plt.figure()
m = Basemap(llcrnrlon=22.5,llcrnrlat=38.5,urcrnrlon=32.,urcrnrlat=43.5,
            rsphere=(6378137.00,6356752.3142),
            resolution='h',projection='merc',
            lat_0=40.,lon_0=20.,lat_ts=20.)

m.drawcoastlines(linewidth=0.2);
m.fillcontinents(color='grey');
m.drawparallels(np.arange(38,44,1),labels=[1,0,0,0]);
m.drawmeridians(np.arange(22,33,2),labels=[0,0,0,1]);

longitude,latitude = m(np.copy(gr.longitude),np.copy(gr.latitude))

im1 = plt.tripcolor(longitude,latitude,gr.element_index-1,MKE
             ,vmin=0,vmax=0.07,shading='gouraud',cmap='needJet2');

cb = m.colorbar(im1,"right", size="5%", pad="4%", extend='max')
cb.set_label('[$m^2/s^2$]',rotation=0,y=1.07,labelpad=-35) 
plt.savefig('paperfigs/uTSS_mean_MKE.png', bbox_inches='tight',format='png',dpi=300)

plt.figure()
m = Basemap(llcrnrlon=22.5,llcrnrlat=38.5,urcrnrlon=32.,urcrnrlat=43.5,
            rsphere=(6378137.00,6356752.3142),
            resolution='h',projection='merc',
            lat_0=40.,lon_0=20.,lat_ts=20.)

m.drawcoastlines(linewidth=0.2);
m.fillcontinents(color='grey');
m.drawparallels(np.arange(38,44,1),labels=[1,0,0,0]);
m.drawmeridians(np.arange(22,33,2),labels=[0,0,0,1]);

longitude,latitude = m(np.copy(gr.longitude),np.copy(gr.latitude))

im1 = plt.tripcolor(longitude,latitude,gr.element_index-1,EKE
             ,vmin=0,vmax=0.07,shading='gouraud',cmap='needJet2');

cb = m.colorbar(im1,"right", size="5%", pad="4%", extend='max')
cb.set_label('[$m^2/s^2$]',rotation=0,y=1.07,labelpad=-35) 
plt.savefig('paperfigs/uTSS_mean_EKE.png', bbox_inches='tight',format='png',dpi=300)


m = Basemap(llcrnrlon=22.5,llcrnrlat=38.5,urcrnrlon=32.,urcrnrlat=43.5,
            rsphere=(6378137.00,6356752.3142),
            resolution='h',projection='merc',
            lat_0=40.,lon_0=20.,lat_ts=20.)
x,y = m(22.7,43.15)

fig, axes = plt.subplots(figsize=(14,6))
ax1 = plt.subplot2grid(shape=(1,2),loc=(0,0), colspan=1)
ax2 = plt.subplot2grid(shape=(1,2),loc=(0,1), colspan=1)
plt.tight_layout()

m.ax = ax1
m.drawcoastlines(linewidth=0.2);
m.fillcontinents(color='grey');
m.drawparallels(np.arange(38,44,1),labels=[0,0,0,0]);
m.drawmeridians(np.arange(22,33,2),labels=[0,0,0,1]);
im1 = ax1.tripcolor(longitude,latitude,gr.element_index-1,MKE
             # ,vmin=0,vmax=0.15,shading='gouraud',cmap='needJet2');
             ,vmin=0,vmax=0.07,shading='gouraud',cmap='needJet2');
ax1.text(x,y,'a)',fontsize=14)
cb = m.colorbar(im1,"left", size="5%", pad="4%", extend='max')
cb.set_label('($m^2/s^2$)',rotation=0,y=1.07,labelpad=-10) 
cb.ax.yaxis.set_ticks_position("left") 

m.ax = ax2
m.drawcoastlines(linewidth=0.2);
m.fillcontinents(color='grey');
m.drawparallels(np.arange(38,44,1),labels=[1,0,0,0]);
m.drawmeridians(np.arange(22,33,2),labels=[0,0,0,1]);
im1 = ax2.tripcolor(longitude,latitude,gr.element_index-1,EKE
             ,vmin=0,vmax=0.07,shading='gouraud',cmap='needJet2');
ax2.text(x,y,'b)',fontsize=14)

cb = m.colorbar(im1,"right", size="5%", pad="4%", extend='max')
cb.set_label('($m^2/s^2$)',rotation=0,y=1.07,labelpad=-35) 
plt.savefig('paperfigs/uTSS_mean_MKE_EKE.png', bbox_inches='tight',format='png',dpi=300)

