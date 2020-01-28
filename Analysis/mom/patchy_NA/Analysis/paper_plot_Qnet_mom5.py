import numpy as np
import sys
import os
# from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt
import scipy.io as io 
import xarray as xr
plt.ion()


root_folder = '/shared/projects/uniklima/globclim/milicak/mom/'

expid = 'om3_core3_ctrl'
gr = xr.open_dataset('/shared/projects/uniklima/globclim/milicak/Analysis/mom/patchy_NA/Analysis/levitus_ewg_temp_salt_cm2m.nc')
# long wave radiation
fnames = root_folder + expid + '/history_63-124years/' + 'lw_heat_00630101.ocean_month.nc'
list = sorted(glob.glob(fnames))
df = xr.open_mfdataset(list)['lw_heat']
tmp = df[516:774,:,:]
LW = tmp.mean('time')
# short wave radiation
fnames = root_folder + expid + '/history_63-124years/' + 'swflx_00630101.ocean_month.nc'
list = sorted(glob.glob(fnames))
df = xr.open_mfdataset(list)['swflx']
tmp = df[516:774,:,:]
SW = tmp.mean('time')
# sensible heat flux
fnames = root_folder + expid + '/history_63-124years/' + 'sens_heat_00630101.ocean_month.nc'
list = sorted(glob.glob(fnames))
df = xr.open_mfdataset(list)['sens_heat']
tmp = df[516:774,:,:]
SH = tmp.mean('time')
# latent heat flux
fnames = root_folder + expid + '/history_63-124years/' + 'evap_heat_00630101.ocean_month.nc'
list = sorted(glob.glob(fnames))
df = xr.open_mfdataset(list)['evap_heat']
tmp = df[516:774,:,:]
LH = tmp.mean('time')
Qnet = SH+LH+LW+SW

expid = 'om3_core3_patchy_full_02'
fnames = root_folder + expid + '/history_63-124years/' + 'lw_heat_00630101.ocean_month.nc'
list = sorted(glob.glob(fnames))
df = xr.open_mfdataset(list)['lw_heat']
tmp = df[516:774,:,:]
LW = tmp.mean('time')
# short wave radiation
fnames = root_folder + expid + '/history_63-124years/' + 'swflx_00630101.ocean_month.nc'
list = sorted(glob.glob(fnames))
df = xr.open_mfdataset(list)['swflx']
tmp = df[516:774,:,:]
SW = tmp.mean('time')
# sensible heat flux
fnames = root_folder + expid + '/history_63-124years/' + 'sens_heat_00630101.ocean_month.nc'
list = sorted(glob.glob(fnames))
df = xr.open_mfdataset(list)['sens_heat']
tmp = df[516:774,:,:]
SH = tmp.mean('time')
# latent heat flux
fnames = root_folder + expid + '/history_63-124years/' + 'evap_heat_00630101.ocean_month.nc'
list = sorted(glob.glob(fnames))
df = xr.open_mfdataset(list)['evap_heat']
tmp = df[516:774,:,:]
LH = tmp.mean('time')
Qnet2 = SH+LH+LW+SW

xloc = -170
yloc = -83
# proj = ccrs.PlateCarree()
disc = plt.get_cmap('Spectral_r',32) 

fig, axes = plt.subplots(figsize=(10,15)) #, subplot_kw=dict(projection=proj))
# plot control
ax1 = plt.subplot(3, 1, 1, projection=ccrs.PlateCarree(central_longitude=-100))
sstbias_plot = plt.pcolormesh(gr.x_T,gr.y_T,Qnet,
               vmin=-200,vmax=200,cmap=disc,transform=ccrs.PlateCarree()
                             ,rasterized=True);
ax1.add_feature(cartopy.feature.LAND,color='grey');
ax1.coastlines();
plt.text(xloc,yloc,'a)',fontsize=18);
# add colorbar
axpos = ax1.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(sstbias_plot, cax=cbar_ax)
cbar.ax.tick_params(labelsize=16)
cbar.set_label('Qnet [W/m^2]', fontsize=16)
ax1.set_yticks([-90,-60,-30,0,30,60,90])
ax1.set_yticklabels([-90,-60,-30,0,30,60,90],fontsize=16) 
ax1.set_ylabel('Lat', fontsize=16)
# plot patchy
ax2 = plt.subplot(3, 1, 2, projection=ccrs.PlateCarree(central_longitude=-100))
sstbias_plot = plt.pcolormesh(gr.x_T,gr.y_T,Qnet2,
               vmin=-200,vmax=200,cmap=disc,transform=ccrs.PlateCarree()
                              ,rasterized=True);
ax2.add_feature(cartopy.feature.LAND,color='grey');
ax2.coastlines();
plt.text(xloc,yloc,'b)',fontsize=18);
# add colorbar
axpos = ax2.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(sstbias_plot, cax=cbar_ax)
cbar.ax.tick_params(labelsize=16)
cbar.set_label('Qnet [W/m^2]', fontsize=16)
ax2.set_yticks([-90,-60,-30,0,30,60,90])
ax2.set_yticklabels([-90,-60,-30,0,30,60,90],fontsize=16) 
ax2.set_ylabel('Lat', fontsize=16)
# plot difference between patchy and control
ax3 = plt.subplot(3, 1, 3, projection=ccrs.PlateCarree(central_longitude=-100))
sstbias_plot = plt.pcolormesh(gr.x_T,gr.y_T,(Qnet2-Qnet),
               vmin=-100,vmax=100,cmap='RdBu_r',transform=ccrs.PlateCarree()
                              ,rasterized=True);
ax3.add_feature(cartopy.feature.LAND,color='grey');
ax3.coastlines();
plt.text(xloc,yloc,'c)',fontsize=18);
# add colorbar
axpos = ax3.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(sstbias_plot, cax=cbar_ax)
cbar.ax.tick_params(labelsize=16)
cbar.set_label('Qnet [W/m^2]', fontsize=16)
ax3.set_xticks(np.array([-180,-120,-60,0,60,120,180]))
ax3.set_yticks([-90,-60,-30,0,30,60,90])
ax3.set_xticklabels(np.array([-180,-120,-60,0,60,120,180])-100,fontsize=16) 
ax3.set_yticklabels([-90,-60,-30,0,30,60,90],fontsize=16) 
ax3.set_xlabel('Lon', fontsize=16)
ax3.set_ylabel('Lat', fontsize=16)

plt.savefig('paperfigs/Qnet_diff_new.png', bbox_inches='tight',format='png',dpi=300)

