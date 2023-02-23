import numpy as np
import glob
from mpl_toolkits.basemap import Basemap
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

ds = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/SSHUGVG/global_ssh_anomaly.nc')

fig, ax1 = plt.subplots(figsize=(8,8))
ax = plt.subplot(1, 1, 1, projection=ccrs.Mercator(central_longitude=0,globe=None))
# ax = plt.axes(projection=ccrs.Mercator(central_longitude=0,globe=None))
# Draw coastlines so we know where we are:
ax.coastlines(resolution='50m', color='black')
# Set the map extent, making sure to specify the correct coordinate system
# for geographical coordinates:
ax.set_extent([-100, 50, -20, 75], crs=ccrs.PlateCarree());
ax.add_feature(cartopy.feature.LAND, facecolor='grey')
im1 = plt.pcolormesh(ds.longitude,ds.latitude,ds['std'],
               vmin=0,vmax=0.35,shading='goaroud',transform=ccrs.PlateCarree()   
                             ,rasterized=True);                            
# add colorbar                                                             
axpos = ax1.get_position()                                                 
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0+0.035,0.03,axpos.height*0.91]);
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

plt.savefig('paperfigs/SSH_std_obs.png', bbox_inches='tight',format='png',dpi=300)

# ax1.set_yticks([-90,-60,-30,0,30,60,90])                                   
# ax1.set_yticklabels([-90,-60,-30,0,30,60,90],fontsize=16)                  
# ax.set_xlabel('Lon', fontsize=16)                                         
# ax.set_ylabel('Lat', fontsize=16)                                         
#

# m = Basemap(projection='merc',llcrnrlat=-20,urcrnrlat=75,
#             llcrnrlon=-100,urcrnrlon=50,lat_ts=20,resolution='i')
# m.drawcoastlines();
# m.fillcontinents(color='grey',lake_color='aqua');
# im1 = m.pcolormesh(gr.geolon,gr.geolat,df.mean_sst,shading='gouraud',
#         vmin=-1,vmax=30,cmap='nice_gfdl',latlon=True);
# cb = m.colorbar(im1,"right", size="5%", pad="2%")
# cb.ax.tick_params(labelsize=14)
#
# parallels = np.arange(25.,50.,5.)
# m.drawparallels(parallels,labels=[True,False,False,False],fontsize=14);
# meridians = np.arange(-90.,-45.,10.)
# m.drawmeridians(meridians,labels=[False,False,False,True],fontsize=14);



