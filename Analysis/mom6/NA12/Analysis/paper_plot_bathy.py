import numpy as np
import glob
from mpl_toolkits.basemap import Basemap
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import xarray as xr

root_folder = '/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/'
gr = xr.open_dataset(root_folder + 'ocean_geometry.nc')
bathy = gr.D*gr.wet



fig, ax1 = plt.subplots(figsize=(8,8))
ax = plt.subplot(1, 1, 1, projection=ccrs.Mercator(central_longitude=0,globe=None))
# ax = plt.axes(projection=ccrs.Mercator(central_longitude=0,globe=None))
# Draw coastlines so we know where we are:
ax.coastlines(resolution='50m', color='black')
# Set the map extent, making sure to specify the correct coordinate system
# for geographical coordinates:
ax.set_extent([-105, 55, -30, 80], crs=ccrs.PlateCarree());
ax.add_feature(cartopy.feature.LAND, facecolor='grey')
im1 = plt.pcolormesh(gr.geolon,gr.geolat,-bathy.where(bathy!=0),
               vmin=-7000,vmax=500,shading='flat',cmap='GMT_ocean',transform=ccrs.PlateCarree()   
                             ,rasterized=True);                            
# add colorbar                                                             
axpos = ax1.get_position()                                                 
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0+0.001,0.03,axpos.height]);
cbar = fig.colorbar(im1, cax=cbar_ax)                             
cbar.ax.tick_params(labelsize=16)                                          
cbar.set_label('m', fontsize=16)                      
lon = np.arange(-100,50,20)
lat = np.array([-20,0,20,40,60,70])
ax.set_xticks(lon, crs=ccrs.PlateCarree())
ax.set_yticks(lat, crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter()
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)

plt.savefig('paperfigs/bathy_MOM6.png', bbox_inches='tight',format='png',dpi=300)
