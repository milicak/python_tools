import numpy as np
import glob
from mpl_toolkits.basemap import Basemap
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

root_folder = '/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/'
gr = xr.open_dataset(root_folder + 'ocean_geometry.nc')

# run this once
# root_folder = '/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/uprime/'
# ls1 = sorted(glob.glob(root_folder+'*uprime*'))
# # last 20 years
# df = xr.open_mfdataset(ls1[8:])
# EKE = df.uprime**2+df.vprime**2
# EKE = EKE.mean('time')
# ds = EKE.to_dataset(name='EKE')
# ds.to_netcdf('/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/uprime/EKE_MOM6.nc')
#

ds = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/uprime/EKE_MOM6.nc')

fig, ax1 = plt.subplots(figsize=(8,8))
ax = plt.subplot(1, 1, 1, projection=ccrs.Mercator(central_longitude=0,globe=None))
# ax = plt.axes(projection=ccrs.Mercator(central_longitude=0,globe=None))
# Draw coastlines so we know where we are:
ax.coastlines(resolution='50m', color='black')
# Set the map extent, making sure to specify the correct coordinate system
# for geographical coordinates:
ax.set_extent([-100, 50, -20, 75], crs=ccrs.PlateCarree());
ax.add_feature(cartopy.feature.LAND, facecolor='grey')
im1 = plt.pcolormesh(gr.geolon,gr.geolat,ds.EKE,
               vmin=0,vmax=0.35,shading='goaroud',cmap='RdYlBu_r',transform=ccrs.PlateCarree()   
                             ,rasterized=True);                            
# add colorbar                                                             
axpos = ax1.get_position()                                                 
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0+0.035,0.03,axpos.height*0.91]);
cbar = fig.colorbar(im1, cax=cbar_ax)                             
cbar.ax.tick_params(labelsize=16)                                          
cbar.set_label(r'm$^2$/s$^2$', fontsize=16)
lon = np.arange(-100,50,20)
lat = np.array([-20,0,20,40,60,70])
ax.set_xticks(lon, crs=ccrs.PlateCarree())
ax.set_yticks(lat, crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter()
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)

plt.savefig('paperfigs/EKE_MOM6.png', bbox_inches='tight',format='png',dpi=300)


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
