import numpy as np
import sys
import os
# from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt
import xarray as xr
plt.ion()

climwoa = '/shared/projects/uniklima/globclim/milicak/Analysis/mom/patchy_NA/Analysis/WOA09_ann_theta_cm2m_extrap.nc'
ds = xr.open_dataset(climwoa,decode_times=False)
dzwoa = np.copy(ds.DEPTH_bnds[:,1]-ds.DEPTH_bnds[:,0]);
dzwoa = dzwoa[:30]

# sstclim = np.transpose(np.copy(ds.tempwoa_mom6[:,:,0]))
tempclim = np.copy(ds.POTENTIAL_TEMP[0,:30,:,:])
nz,ny,nx = tempclim.shape 
dzwoa3d = np.tile(np.copy(dzwoa),(nx,ny,1))
dzwoa3d = np.transpose(dzwoa3d,[2,1,0])

root_folder = '/shared/projects/uniklima/globclim/milicak/mom/'

expid = 'om3_core3_ctrl'
gr = xr.open_dataset('/shared/projects/uniklima/globclim/milicak/Analysis/mom/patchy_NA/Analysis/levitus_ewg_temp_salt_cm2m.nc')
fnames = root_folder + expid + '/history_63-124years/' + 'temp_00630101.ocean_month.nc'
list = sorted(glob.glob(fnames))
df = xr.open_mfdataset(list)['temp']
data = df[348:744,:30,:,:]
data = data.mean('time')

expid = 'om3_core3_patchy_full_02'
fnames = root_folder + expid + '/history_63-124years/' + 'temp_00630101.ocean_month.nc'
list = sorted(glob.glob(fnames))
df2 = xr.open_mfdataset(list)['temp']
data2 = df2[348:744,:30,:,:]
data2 = data2.mean('time')

dnm = np.copy(data.data)
dzwoa3d[np.isnan(dnm)] = np.nan
tempclim[np.isnan(dnm)] = np.nan

dnm = tempclim*dzwoa3d                
dnm = np.nansum(dnm,axis=0)       
dz3dsum = np.nansum(dzwoa3d,axis=0)  
tempclim = dnm/dz3dsum            

dnm = np.copy(data.data)*dzwoa3d    
dnm = np.nansum(dnm,axis=0)       
temp1 = dnm/dz3dsum               
                                  
dnm = np.copy(data2.data)*dzwoa3d    
dnm = np.nansum(dnm,axis=0)       
temp2 = dnm/dz3dsum               


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
sstbias_plot = plt.pcolormesh(gr.x_T,gr.y_T,temp1-tempclim,
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
sstbias_plot = plt.pcolormesh(gr.x_T,gr.y_T,temp2-temp1,
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
cbar.ax.tick_params(labelsize=16)
cbar.set_label('Temperature [$\circ$C]', fontsize=16)
ax3.set_yticks([20,40,60,80])                    
ax3.set_yticklabels([20,40,60,80],fontsize=16)   
ax3.set_xticks([-120,-90,-60,-30,0])                            
ax3.set_xticklabels(np.array([-120,-90,-60,-30,0]),fontsize=16) 
ax3.set_xlabel('Lon', fontsize=16)
ax3.set_ylabel('Lat', fontsize=16)

plt.savefig('paperfigs/temp_500m_bias_diff_NA_mom5.png', bbox_inches='tight',format='png',dpi=300)

