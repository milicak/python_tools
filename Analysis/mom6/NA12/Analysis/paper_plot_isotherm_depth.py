import numpy as np
import glob
from mpl_toolkits.basemap import Basemap
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

root_folder = '/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/'
gr = xr.open_dataset(root_folder + 'ocean_geometry.nc')
df = xr.open_dataset('Isotherm_depth_MOM6.nc')
ds = xr.open_dataset('Isotherm_depth_climatology.nc', decode_times=False)


mom6_iso_therm = df.isotherm_depth.mean('time')

fig, ax1 = plt.subplots(figsize=(8,8))
ax = plt.subplot(1, 1, 1, projection=ccrs.Mercator(central_longitude=0,globe=None))
ax.add_feature(cartopy.feature.LAND,zorder=2, facecolor='grey')
ax.coastlines(resolution='50m', zorder=3, color='black')
ax.set_extent([-100, 30, -20, 45], crs=ccrs.PlateCarree());
ax.add_feature(cartopy.feature.LAND, facecolor='grey')
# im1 = plt.pcolormesh(gr.geolon,gr.geolat,mom6_iso_therm,
im1 = plt.pcolormesh(gr.geolon,gr.geolat,mom6_iso_therm.where(mom6_iso_therm!=2.5),
        vmin=0,vmax=150,shading='flat',cmap='RdYlBu_r',transform=ccrs.PlateCarree()   
                             ,rasterized=True);                            
# add colorbar                                                             
axpos = ax1.get_position()                                                 
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0+0.056,0.03,axpos.height*0.852]);
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
plt.savefig('paperfigs/mean_isotherm_depth_MOM6.png', bbox_inches='tight',format='png',dpi=300)

fig, ax1 = plt.subplots(figsize=(8,8))
ax = plt.subplot(1, 1, 1, projection=ccrs.Mercator(central_longitude=0,globe=None))
# Draw coastlines so we know where we are:
# Set the map extent, making sure to specify the correct coordinate system
# for geographical coordinates:
ax.set_extent([-100, 50, -20, 75], crs=ccrs.PlateCarree());
ax.add_feature(cartopy.feature.LAND,zorder=2, facecolor='grey')
ax.coastlines(resolution='50m', zorder=3, color='black')
# im1 = plt.pcolormesh(ds.lon,ds.lat,ds.isotherm_depth[0,:,:],
im1 = plt.pcolormesh(ds.lon,ds.lat,ds.isotherm_depth[0,:,:].where(ds.isotherm_depth[0,:,:]!=0),
        vmin=0,vmax=150,shading='flat',cmap='RdYlBu_r',transform=ccrs.PlateCarree()   
                             ,rasterized=True);                            
# add colorbar                                                             
axpos = ax1.get_position()                                                 
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0+0.056,0.03,axpos.height*0.852]);
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
plt.savefig('paperfigs/mean_isotherm_depth_obs.png', bbox_inches='tight',format='png',dpi=300)

