import numpy as np
import sys
import os
# from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt
import xarray as xr
plt.ion()


root_folder = '/archive/milicak/MOM6-examples/Projects/patchy_NA/'

expid = 'work_patchy'
gr = xr.open_dataset(root_folder+expid+'/ocean_geometry.nc')
fnames = root_folder + expid + '/OUT/' + '19050101.ocean_month.nc'
list = sorted(glob.glob(fnames))
fnames = root_folder + expid + '/OUT/' + '19100101.ocean_month.nc'
list.extend(sorted(glob.glob(fnames)))
# df = xr.open_mfdataset(list)['MEKE_deltarho_TW']
df = xr.open_mfdataset(list)['MEKE_deltarho_eos']
# data = df[:,:,:]
data = df[:,0,:,:]
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
# ax1 = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=-100.))
ax1 = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
sstbias_plot = plt.pcolormesh(gr.geolon,gr.geolat,3*ma.masked_equal(data,0),
               vmin=-.5,vmax=0,shading='flat',
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
# ax1.set_xticklabels(np.array([-180,-120,-60,0,60,120,180])-100,fontsize=16)
ax1.set_xticklabels(np.array([-180,-120,-60,0,60,120,180]),fontsize=16)
ax1.set_xlabel('Lon', fontsize=16)
ax1.set_ylabel('Lat', fontsize=16)

plt.savefig('paperfigs/deltarho.png', bbox_inches='tight',format='png',dpi=300)

