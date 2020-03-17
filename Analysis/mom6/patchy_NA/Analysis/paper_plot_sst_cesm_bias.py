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

archive/milicak/cesm_simulations/GMOM_NApatchy_eos/run/GMOM_NApatchy_eos.mom6.hm_0033_04.nc


fnames = root_folder + expid + '/run/' + 'GMOM_NApatchy_eos.mom6.hm_003*.nc'
list = sorted(glob.glob(fnames))
df = xr.open_mfdataset(list)['tos']
data = df[:,:,:]
data = data.mean('time')

# ax = plt.axes(projection=ccrs.Robinson());
# plt.pcolormesh(gr.geolon,gr.geolat,df.data[6,0,:,:],transform=ccrs.Robinson());plt.colorbar();
# ax = plt.axes(projection=ccrs.PlateCarree())
# plt.pcolormesh(gr.geolon,gr.geolat,data,transform=ccrs.PlateCarree());plt.colorbar();

proj = ccrs.PlateCarree()

fig, axes = plt.subplots(figsize=(10,15)) #, subplot_kw=dict(projection=proj))
# plot control
ax1 = plt.subplot(3, 1, 1, projection=ccrs.PlateCarree())
sstbias_plot = plt.pcolormesh(gr.TLONG,gr.TLAT,data-sstclim,
               vmin=-5,vmax=5,cmap='RdBu_r',transform=ccrs.PlateCarree()
                             ,rasterized=True);
ax1.add_feature(cartopy.feature.LAND,color='grey');
ax1.coastlines();
# add colorbar
axpos = ax1.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(sstbias_plot, cax=cbar_ax)
cbar.ax.tick_params(labelsize=12)
cbar.set_label('Temperature ($\circ$C)', fontsize=12)
# plot patchy
ax2 = plt.subplot(3, 1, 2, projection=ccrs.PlateCarree())
sstbias_plot = plt.pcolormesh(gr.TLONG,gr.TLAT,data2-sstclim,
               vmin=-5,vmax=5,cmap='RdBu_r',transform=ccrs.PlateCarree()
                              ,rasterized=True);
ax2.add_feature(cartopy.feature.LAND,color='grey');
ax2.coastlines();
# add colorbar
axpos = ax2.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(sstbias_plot, cax=cbar_ax)
cbar.ax.tick_params(labelsize=12)
cbar.set_label('Temperature ($\circ$C)', fontsize=12)
# plot difference between patchy and control
ax3 = plt.subplot(3, 1, 3, projection=ccrs.PlateCarree())
sstbias_plot = plt.pcolormesh(gr.TLON,gr.TLAT,data2-data,
               vmin=-2.5,vmax=2.5,cmap='RdBu_r',transform=ccrs.PlateCarree()
                              ,rasterized=True);
ax3.add_feature(cartopy.feature.LAND,color='grey');
ax3.coastlines();
# add colorbar
axpos = ax3.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(sstbias_plot, cax=cbar_ax)
cbar.ax.tick_params(labelsize=12)
cbar.set_label('Temperature ($\circ$C)', fontsize=12)

plt.savefig('paperfigs/sst_cesm_bias_diff.png', bbox_inches='tight',format='png',dpi=300)

