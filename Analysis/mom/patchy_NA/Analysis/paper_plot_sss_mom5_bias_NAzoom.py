import numpy as np
import sys
import os
# from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt
import xarray as xr
plt.ion()

climwoa = '/shared/projects/uniklima/globclim/milicak/Analysis/mom/patchy_NA/Analysis/WOA09_ann_salinity_cm2m_extrap.nc'
ds = xr.open_dataset(climwoa,decode_times=False)
# sstclim = np.transpose(np.copy(ds.tempwoa_mom6[:,:,0]))
sssclim = np.copy(ds.S_AN[0,0,:,:])

root_folder = '/shared/projects/uniklima/globclim/milicak/mom/'
xoffset = -100

expid = 'om3_core3_ctrl'
gr = xr.open_dataset('/shared/projects/uniklima/globclim/milicak/Analysis/mom/patchy_NA/Analysis/levitus_ewg_temp_salt_cm2m.nc')
fnames = root_folder + expid + '/history_63-124years/' + 'salt_00630101.ocean_month.nc'
list = sorted(glob.glob(fnames))
df = xr.open_mfdataset(list)['salt']
data = df[348:744,0,:,:]
data = data.mean('time')

expid = 'om3_core3_patchy_full_02'
fnames = root_folder + expid + '/history_63-124years/' + 'salt_00630101.ocean_month.nc'
list = sorted(glob.glob(fnames))
df2 = xr.open_mfdataset(list)['salt']
data2 = df2[348:744,0,:,:]
data2 = data2.mean('time')

# proj = ccrs.PlateCarree()
xloc = -170
yloc = -83

fig, axes = plt.subplots(figsize=(10,10)) #, subplot_kw=dict(projection=proj))
# fig.tight_layout()
# plot control
# ax1 = plt.subplot(2, 1, 1, projection=ccrs.PlateCarree(central_longitude=-100))
ax1 = plt.subplot(2, 1, 1, projection=ccrs.PlateCarree())
sstbias_plot = plt.pcolormesh(gr.x_T,gr.y_T,data-sssclim,
               vmin=-3,vmax=3,cmap='RdBu_r',transform=ccrs.PlateCarree()
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
cbar.set_label('Salinity [psu]', fontsize=16)
ax1.set_yticks([20,40,60,80])                    
ax1.set_yticklabels([20,40,60,80],fontsize=16)   
ax1.set_ylabel('Lat', fontsize=16)
# plot patchy
ax3 = plt.subplot(2, 1, 2,projection=ccrs.PlateCarree())
sstbias_plot = plt.pcolormesh(gr.x_T,gr.y_T,data2-data,
               vmin=-1,vmax=1,cmap='RdBu_r',transform=ccrs.PlateCarree()
                              ,rasterized=True);
ax3.add_feature(cartopy.feature.LAND,color='grey');
ax3.coastlines();
ax3.set_extent([-120, 15, 10, 80], ccrs.PlateCarree())
axpos = ax3.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(sstbias_plot, cax=cbar_ax)
cbar.ax.tick_params(labelsize=16)
cbar.set_label('Salinity [psu]', fontsize=16)
ax3.set_yticks([20,40,60,80])                    
ax3.set_yticklabels([20,40,60,80],fontsize=16)   
ax3.set_xticks([-120,-90,-60,-30,0])                            
ax3.set_xticklabels(np.array([-120,-90,-60,-30,0]),fontsize=16) 
ax3.set_xlabel('Lon', fontsize=16)
ax3.set_ylabel('Lat', fontsize=16)

plt.savefig('paperfigs/sss_bias_diff_NA_mom5.png', bbox_inches='tight',format='png',dpi=300)

