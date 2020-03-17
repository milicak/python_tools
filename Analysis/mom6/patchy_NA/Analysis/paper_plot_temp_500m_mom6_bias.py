import numpy as np
import sys
import os
from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt
import xarray as xr
plt.ion()

climwoa = 'tempwoa_mom6_0_5deg.nc'
ds = xr.open_dataset(climwoa)
sstclim = np.transpose(np.copy(ds.tempwoa_mom6[:,:,:14]))
nz,ny,nx = sstclim.shape

root_folder = '/archive/milicak/MOM6-examples/Projects/patchy_NA/'

expid = 'work_ctrl'
gr = xr.open_dataset(root_folder+expid+'/ocean_geometry.nc')
fnames = root_folder + expid + '/OUT/' + '19050101.ocean_month_z.nc'
list = sorted(glob.glob(fnames))
fnames = root_folder + expid + '/OUT/' + '19100101.ocean_month_z.nc'
list.extend(sorted(glob.glob(fnames)))
df = xr.open_mfdataset(list)['thetao']
data = df[:,:14,:,:]
data = data.mean('time')
zw = np.copy(xr.open_dataset(fnames)['z_i'])
dz = zw[1:]-zw[0:-1]
dz3d = np.tile(dz[:14],[nx, ny, 1])
dz3d = np.transpose(dz3d,[2,1,0])
dz3d[np.isnan(sstclim)] = np.nan

dnm = sstclim*dz3d
dnm = np.nansum(dnm,axis=0)
dz3dsum = np.nansum(dz3d,axis=0)
tempclim = dnm/dz3dsum

expid = 'work_patchy'
fnames = root_folder + expid + '/OUT/' + '19050101.ocean_month_z.nc'
list = sorted(glob.glob(fnames))
fnames = root_folder + expid + '/OUT/' + '19100101.ocean_month_z.nc'
list.extend(sorted(glob.glob(fnames)))
df2 = xr.open_mfdataset(list)['thetao']
data2 = df2[:,:14,:,:]
data2 = data2.mean('time')

dnm = np.copy(data.data)*dz3d
dnm = np.nansum(dnm,axis=0)
temp1 = dnm/dz3dsum

dnm = np.copy(data2.data)*dz3d
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
sstbias_plot = plt.pcolormesh(gr.geolon,gr.geolat,temp1-tempclim,
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
# ax1.set_yticks([-90,-60,-30,0,30,60,90])
# ax1.set_yticklabels([-90,-60,-30,0,30,60,90],fontsize=16)
ax1.set_yticks([20,40,60,80])
ax1.set_yticklabels([20,40,60,80],fontsize=16)
ax1.set_ylabel('Lat', fontsize=16)
# plot patchy
ax3 = plt.subplot(2, 1, 2, projection=ccrs.PlateCarree())
sstbias_plot = plt.pcolormesh(gr.geolon,gr.geolat,temp2-temp1,
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
ax3.set_xticks([-120,-90,-60,-30,0])
ax3.set_xticklabels(np.array([-120,-90,-60,-30,0]),fontsize=16)
ax3.set_yticks([20,40,60,80])
ax3.set_yticklabels([20,40,60,80],fontsize=16)
ax3.set_xlabel('Lon', fontsize=16)
ax3.set_ylabel('Lat', fontsize=16)

plt.savefig('paperfigs/temp500m_bias_diff_NA_mom6.png', bbox_inches='tight',format='png',dpi=300)

