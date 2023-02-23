import numpy as np
import glob
from mpl_toolkits.basemap import Basemap
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

root_folder = '/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/'
gr = xr.open_dataset(root_folder + 'ocean_geometry.nc')
df = xr.open_dataset('GS_separation_mom6.nc')
ds = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/sst_woa_NA12.nc',decode_times=False)

fig, ax1 = plt.subplots(figsize=(8,8))
ax = plt.subplot(1, 1, 1, projection=ccrs.Mercator(central_longitude=0,globe=None))
# ax = plt.axes(projection=ccrs.Mercator(central_longitude=0,globe=None))
# Draw coastlines so we know where we are:
ax.coastlines(resolution='50m', color='black')
# Set the map extent, making sure to specify the correct coordinate system
# for geographical coordinates:
ax.set_extent([-100, 50, -20, 75], crs=ccrs.PlateCarree());
ax.add_feature(cartopy.feature.LAND, facecolor='grey')
im1 = plt.pcolormesh(gr.geolon,gr.geolat,df.mean_sst,
        vmin=-1,vmax=32,shading='goaroud',cmap='nice_gfdl',transform=ccrs.PlateCarree()   
                             ,rasterized=True);                            
# add colorbar                                                             
axpos = ax1.get_position()                                                 
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0+0.035,0.03,axpos.height*0.91]);
cbar = fig.colorbar(im1, cax=cbar_ax)                             
cbar.ax.tick_params(labelsize=16)                                          
# cbar.set_label('m/s', fontsize=16)                      
lon = np.arange(-100,50,20)
lat = np.array([-20,0,20,40,60,70])
ax.set_xticks(lon, crs=ccrs.PlateCarree())
ax.set_yticks(lat, crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter()
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
plt.savefig('paperfigs/mean_sst_MOM6.png', bbox_inches='tight',format='png',dpi=300)

fig, ax1 = plt.subplots(figsize=(8,8))
ax = plt.subplot(1, 1, 1, projection=ccrs.Mercator(central_longitude=0,globe=None))
# Draw coastlines so we know where we are:
ax.coastlines(resolution='50m', color='black')
# Set the map extent, making sure to specify the correct coordinate system
# for geographical coordinates:
ax.set_extent([-100, 50, -20, 75], crs=ccrs.PlateCarree());
ax.add_feature(cartopy.feature.LAND, facecolor='grey')
im1 = plt.pcolormesh(gr.geolon,gr.geolat,df.mean_sst-np.copy(ds.sst_woa),
        vmin=-4,vmax=4,shading='goaroud',cmap='RdBu_r',transform=ccrs.PlateCarree()   
                             ,rasterized=True);                            
# add colorbar                                                             
axpos = ax1.get_position()                                                 
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0+0.035,0.03,axpos.height*0.91]);
cbar = fig.colorbar(im1, cax=cbar_ax)                             
cbar.ax.tick_params(labelsize=16)                                          
# cbar.set_label('m/s', fontsize=16)                      
lon = np.arange(-100,50,20)
lat = np.array([-20,0,20,40,60,70])
ax.set_xticks(lon, crs=ccrs.PlateCarree())
ax.set_yticks(lat, crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter()
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
plt.savefig('paperfigs/mean_sst_bias_MOM6.png', bbox_inches='tight',format='png',dpi=300)

