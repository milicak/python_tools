import numpy as np
import sys
import os
from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt
import xarray as xr
plt.ion()

climwoa = '/archive/milicak/dataset/CESM/inputdata/ocn/mom/tx0.66v1/PHC2_TEMP_tx0.66v1_34lev_ann_avg.nc'
ds = xr.open_dataset(climwoa)
sstclim = np.copy(ds.TEMP[0,:,:])

root_folder = '/archive/milicak/cesm_simulations/'

expid = 'GMOM_NApatchy_eos'
gr = xr.open_dataset(climwoa)

fnames = root_folder + expid + '/run/' + 'GMOM_NApatchy_eos.mom6.hm_003*.nc'
list = sorted(glob.glob(fnames))
list = list[:108]
df = xr.open_mfdataset(list)['tos']
data = df[:,:,:]
data = data.mean('time')

expid = 'GMOM_NApatchy_eos3'
fnames = root_folder + expid + '/run/' + 'GMOM_NApatchy_eos3.mom6.hm_003*.nc'
list = sorted(glob.glob(fnames))
list = list[:108]
df2 = xr.open_mfdataset(list)['tos']
data2 = df2[:,:,:]
data2 = data2.mean('time')

# ax = plt.axes(projection=ccrs.Robinson());
# plt.pcolormesh(gr.geolon,gr.geolat,df.data[6,0,:,:],transform=ccrs.Robinson());plt.colorbar();
# ax = plt.axes(projection=ccrs.PlateCarree())
# plt.pcolormesh(gr.geolon,gr.geolat,data,transform=ccrs.PlateCarree());plt.colorbar();

xloc = -170
yloc = -83
proj = ccrs.PlateCarree()

fig, axes = plt.subplots(figsize=(10,10)) #, subplot_kw=dict(projection=proj))
# plot control
ax1 = plt.subplot(2, 1, 1, projection=ccrs.PlateCarree())
sstbias_plot = plt.pcolormesh(gr.TLONG,gr.TLAT,data-sstclim,
               vmin=-5,vmax=5,cmap='RdBu_r',transform=ccrs.PlateCarree()
                             ,rasterized=True);
ax1.add_feature(cartopy.feature.LAND,color='grey');
ax1.coastlines();
ax1.set_extent([-120, 15, 10, 80], ccrs.PlateCarree())
# plt.text(xloc,yloc,'a)',fontsize=18);
# add colorbar
axpos = ax1.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(sstbias_plot, cax=cbar_ax)
cbar.ax.tick_params(labelsize=16)
cbar.set_label('Temperature [$\circ$C]', fontsize=16)
ax1.set_yticks([20,40,60,80])
ax1.set_yticklabels([20,40,60,80],fontsize=16)
ax1.set_ylabel('Lat', fontsize=16)
# plot patchy
ax3 = plt.subplot(2, 1, 2, projection=ccrs.PlateCarree())
sstbias_plot = plt.pcolormesh(gr.TLONG,gr.TLAT,data2-data,
               vmin=-2.5,vmax=2.5,cmap='RdBu_r',transform=ccrs.PlateCarree()
                              ,rasterized=True);
ax3.add_feature(cartopy.feature.LAND,color='grey');
ax3.coastlines();
ax3.set_extent([-120, 15, 10, 80], ccrs.PlateCarree())
# plt.text(xloc,yloc,'c)',fontsize=18);
# add colorbar
axpos = ax3.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(sstbias_plot, cax=cbar_ax)
cbar.ax.tick_params(labelsize=12)
cbar.set_label('Temperature [$\circ$C]', fontsize=12)
ax3.set_xticks([-120,-90,-60,-30,0])
ax3.set_xticklabels(np.array([-120,-90,-60,-30,0]),fontsize=16)
ax3.set_yticks([20,40,60,80])
ax3.set_yticklabels([20,40,60,80],fontsize=16)
ax3.set_xlabel('Lon', fontsize=16)
ax3.set_ylabel('Lat', fontsize=16)

plt.savefig('paperfigs/sst_cesm_bias_diff_NA.png', bbox_inches='tight',format='png',dpi=300)

