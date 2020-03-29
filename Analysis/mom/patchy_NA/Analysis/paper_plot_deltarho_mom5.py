import numpy as np
import sys
import os
# from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt
import xarray as xr
plt.ion()


root_folder = '/shared/projects/uniklima/globclim/milicak/mom/'

gr = xr.open_dataset('/shared/projects/uniklima/globclim/milicak/Analysis/mom/patchy_NA/Analysis/levitus_ewg_temp_salt_cm2m.nc')

expid = 'om3_core3_patchy_full_02'
fnames = root_folder + expid + '/history_63-124years/' + 'delta_rho_anom_00630101.ocean_month.nc'
# fnames = root_folder + expid + '/history_63-124years/' + 'patchy_delta_rho_tw_00630101.ocean_month.nc'
list = sorted(glob.glob(fnames))
df = xr.open_mfdataset(list)['delta_rho_anom']
# df = xr.open_mfdataset(list)['patchy_delta_rho_tw']
data = df[348:744,0,:,:]
# data = df[348:744,:,:]
data = data.mean('time')

# ax = plt.axes(projection=ccrs.Robinson());
# plt.pcolormesh(gr.geolon,gr.geolat,df.data[6,0,:,:],transform=ccrs.Robinson());plt.colorbar();
# ax = plt.axes(projection=ccrs.PlateCarree())
# plt.pcolormesh(gr.geolon,gr.geolat,data,transform=ccrs.PlateCarree());plt.colorbar();

xloc = -170 
yloc = -83  

proj = ccrs.PlateCarree()
disc = plt.get_cmap('RdBu',32);
# disc = plt.cm.get_cmap('RdBu_r',32);

fig, axes = plt.subplots(figsize=(10,5)) #, subplot_kw=dict(projection=proj))
# ax1 = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax1 = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=-100))
sstbias_plot = plt.pcolormesh(gr.x_T,gr.y_T,ma.masked_equal(data,0),
               vmin=-1,vmax=0,shading='flat',
               # vmin=0,vmax=2,shading='flat',
               cmap=disc,transform=ccrs.PlateCarree()
                             ,rasterized=True);
ax1.add_feature(cartopy.feature.LAND,color='grey');
ax1.coastlines();
# plt.text(xloc,yloc,'a)',fontsize=18); 
# add colorbar
axpos = ax1.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(sstbias_plot, cax=cbar_ax)
cbar.ax.tick_params(labelsize=16)
# cbar.set_label('$\Delta \rho$ [kg/m-3]', fontsize=16)
ax1.set_yticks([-90,-60,-30,0,30,60,90])
ax1.set_yticklabels([-90,-60,-30,0,30,60,90],fontsize=16) 
ax1.set_xticks([-180,-120,-60,0,60,120,180])
# ax1.set_xticklabels(np.array([-180,-120,-60,0,60,120,180]),fontsize=16) 
ax1.set_xticklabels(np.array([-180,-120,-60,0,60,120,180])-100,fontsize=16) 
ax1.set_xlabel('Lon', fontsize=16)
ax1.set_ylabel('Lat', fontsize=16)

plt.savefig('paperfigs/deltarho.png', bbox_inches='tight',format='png',dpi=300)

