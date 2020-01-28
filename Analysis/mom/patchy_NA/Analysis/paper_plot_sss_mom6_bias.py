import numpy as np
import sys
import os
from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt
import xarray as xr
plt.ion()

climwoa = 'saltwoa_mom6_0_5deg.nc'
ds = xr.open_dataset(climwoa)
sssclim = np.transpose(np.copy(ds.saltwoa_mom6[:,:,0]))

root_folder = '/archive/milicak/MOM6-examples/Projects/patchy_NA/'

expid = 'work_ctrl'
gr = xr.open_dataset(root_folder+expid+'/ocean_geometry.nc')
fnames = root_folder + expid + '/OUT/' + '19050101.ocean_month_z.nc'
list = sorted(glob.glob(fnames))
fnames = root_folder + expid + '/OUT/' + '19100101.ocean_month_z.nc'
list.extend(sorted(glob.glob(fnames)))
df = xr.open_mfdataset(list)['so']
data = df[:,0,:,:]
data = data.mean('time')

expid = 'work_patchy'
fnames = root_folder + expid + '/OUT/' + '19050101.ocean_month_z.nc'
list = sorted(glob.glob(fnames))
fnames = root_folder + expid + '/OUT/' + '19100101.ocean_month_z.nc'
list.extend(sorted(glob.glob(fnames)))
df2 = xr.open_mfdataset(list)['so']
data2 = df2[:,0,:,:]
data2 = data2.mean('time')

# ax = plt.axes(projection=ccrs.Robinson());
# plt.pcolormesh(gr.geolon,gr.geolat,df.data[6,0,:,:],transform=ccrs.Robinson());plt.colorbar();
# ax = plt.axes(projection=ccrs.PlateCarree())
# plt.pcolormesh(gr.geolon,gr.geolat,data,transform=ccrs.PlateCarree());plt.colorbar();

proj = ccrs.PlateCarree()

fig, axes = plt.subplots(figsize=(10,15)) #, subplot_kw=dict(projection=proj))
# plot control
ax1 = plt.subplot(3, 1, 1, projection=ccrs.PlateCarree())
sstbias_plot = plt.pcolormesh(gr.geolon,gr.geolat,data-sssclim,
               vmin=-3,vmax=3,cmap='RdBu_r',transform=ccrs.PlateCarree()
                             ,rasterized=True);
ax1.add_feature(cartopy.feature.LAND,color='grey');
ax1.coastlines();
# add colorbar
axpos = ax1.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(sstbias_plot, cax=cbar_ax)
cbar.ax.tick_params(labelsize=12)
cbar.set_label('Salinity (psu)', fontsize=12)
# plot patchy
ax2 = plt.subplot(3, 1, 2, projection=ccrs.PlateCarree())
sstbias_plot = plt.pcolormesh(gr.geolon,gr.geolat,data2-sssclim,
               vmin=-3,vmax=3,cmap='RdBu_r',transform=ccrs.PlateCarree()
                              ,rasterized=True);
ax2.add_feature(cartopy.feature.LAND,color='grey');
ax2.coastlines();
# add colorbar
axpos = ax2.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(sstbias_plot, cax=cbar_ax)
cbar.ax.tick_params(labelsize=12)
cbar.set_label('Salinity (psu)', fontsize=12)
# plot difference between patchy and control
ax3 = plt.subplot(3, 1, 3, projection=ccrs.PlateCarree())
sstbias_plot = plt.pcolormesh(gr.geolon,gr.geolat,data2-data,
               vmin=-1,vmax=1,cmap='RdBu_r',transform=ccrs.PlateCarree()
                              ,rasterized=True);
ax3.add_feature(cartopy.feature.LAND,color='grey');
ax3.coastlines();
# add colorbar
axpos = ax3.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(sstbias_plot, cax=cbar_ax)
cbar.ax.tick_params(labelsize=12)
cbar.set_label('Salinity (psu)', fontsize=12)

plt.savefig('paperfigs/sss_bias_diff.png', bbox_inches='tight',format='png',dpi=300)

