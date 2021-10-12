import numpy as np
import matplotlib.colors as colors  
from mpl_toolkits.basemap import Basemap
import xarray as xr
import cmocean


root_folder = '/work/opa/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/monthly_clim/'

fname = root_folder + 'uTSS_TS_monthly_clim.nc'

df = xr.open_dataset(fname)

ds = df.salinity.mean('month')
dt = df.temperature.mean('month')
divnorm = colors.DivergingNorm(vmin=17, vcenter=22, vmax=38) 

# salinity
plt.figure()
m = Basemap(llcrnrlon=22.5,llcrnrlat=38.5,urcrnrlon=32.,urcrnrlat=43.5,
            rsphere=(6378137.00,6356752.3142),
            resolution='h',projection='merc',
            lat_0=40.,lon_0=20.,lat_ts=20.)

m.drawcoastlines(linewidth=0.2);
m.fillcontinents(color='grey');
m.drawparallels(np.arange(38,44,1),labels=[1,0,0,0]);
m.drawmeridians(np.arange(22,33,2),labels=[0,0,0,1]);

longitude,latitude = m(np.copy(df.longitude[0,:]),np.copy(df.latitude[0,:]))

im1 = plt.tripcolor(longitude,latitude,df.element_index[-1,:,:]-1,ds[:,0]
             ,shading='gouraud',cmap='needJet2',norm=divnorm);

cb = m.colorbar(im1,"right", size="5%", pad="4%", extend='both')
cb.set_label('psu',rotation=0,y=1.07,labelpad=-36)
plt.savefig('paperfigs/uTSS_mean_salinity.png', bbox_inches='tight',format='png',dpi=300)


plt.figure()
m = Basemap(llcrnrlon=22.5,llcrnrlat=38.5,urcrnrlon=32.,urcrnrlat=43.5,
            rsphere=(6378137.00,6356752.3142),
            resolution='h',projection='merc',
            lat_0=40.,lon_0=20.,lat_ts=20.)

m.drawcoastlines(linewidth=0.2);
m.fillcontinents(color='grey');
m.drawparallels(np.arange(38,44,1),labels=[1,0,0,0]);
m.drawmeridians(np.arange(22,33,2),labels=[0,0,0,1]);

longitude,latitude = m(np.copy(df.longitude[0,:]),np.copy(df.latitude[0,:]))

im1 = plt.tripcolor(longitude,latitude,df.element_index[-1,:,:]-1,dt[:,0]
             ,shading='gouraud',cmap='needJet2',vmin=14,vmax=21);

cb = m.colorbar(im1,"right", size="5%", pad="4%", extend='both')
cb.set_label('$^\circ$C',rotation=0,y=1.07,labelpad=-30)
plt.savefig('paperfigs/uTSS_mean_temperature.png', bbox_inches='tight',format='png',dpi=300)

