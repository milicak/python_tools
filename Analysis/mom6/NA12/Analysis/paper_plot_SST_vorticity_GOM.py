import numpy as np
import glob
from mpl_toolkits.basemap import Basemap
import cartopy
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from xgcm import Grid

def coriolis(lat):
    """Compute the Coriolis parameter for the given latitude:
    ``f = 2*omega*sin(lat)``, where omega is the angular velocity
    of the Earth.

    Parameters
    ----------
    lat : array
      Latitude [degrees].
    """
    omega   = 7.2921159e-05  # angular velocity of the Earth [rad/s]
    return 2*omega*np.sin(lat/360.*2*np.pi)


deg_rad=np.pi/180.

gr = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/ocean_geometry.nc')
gr = gr.rename({'lonh':'xh','lath':'yh','lonq':'xq','latq':'yq'})
latb = np.copy(gr.geolatb)
fcor = coriolis(latb)

ls1 = sorted(glob.glob('/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/ocean_daily/2014*daily*'))
df = xr.open_dataset(ls1[0])
df = df.merge(gr)
# 2D grid
grid = Grid(df, coords={'X': {'center': 'xh', 'outer': 'xq'},
                        'Y': {'center': 'yh', 'outer': 'yq'},
                         }, periodic=[])

vorticity = ( - grid.diff(df.ssu * df.dxCu, 'Y', boundary='fill')
              + grid.diff(df.ssv * df.dyCv, 'X', boundary='fill') ) / df.Aq

fig, ax1 = plt.subplots(figsize=(8,16))                                              
ax = plt.subplot(2, 1, 1, projection=ccrs.Mercator(central_longitude=0,globe=None))
# ax = plt.axes(projection=ccrs.Mercator(central_longitude=0,globe=None))
# Draw coastlines so we know where we are:
ax.coastlines(resolution='50m', color='black')
# Set the map extent, making sure to specify the correct coordinate system
# for geographical coordinates:
ax.set_extent([-100, -30, 10, 55], crs=ccrs.PlateCarree());
ax.add_feature(cartopy.feature.LAND, facecolor='grey')
im1 = plt.pcolormesh(gr.geolon,gr.geolat,df.tos[180,:,:],
        vmin=2,vmax=34,shading='auto',cmap='nice_gfdl',
        transform=ccrs.PlateCarree(),rasterized=True);                            
# add colorbar                                                             
lon = np.arange(-100,-20,20)
lat = np.array([10,20,30,40,50])
ax.set_xticks(lon, crs=ccrs.PlateCarree())
ax.set_yticks(lat, crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter()
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)

axpos = ax.get_position()                                                 
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0+0.005,0.03,axpos.height*0.98]);
cbar = fig.colorbar(im1, cax=cbar_ax)                             
cbar.ax.tick_params(labelsize=16)                                          
cbar.set_label(r'$^{\circ}$C', fontsize=16)                      

ax = plt.subplot(2, 1, 2, projection=ccrs.Mercator(central_longitude=0,globe=None))  
ax.coastlines(resolution='50m', color='black')
ax.set_extent([-100, -30, 10, 55], crs=ccrs.PlateCarree());
ax.add_feature(cartopy.feature.LAND, facecolor='grey')
im1 = plt.pcolormesh(gr.geolonb,gr.geolatb,vorticity[180,:,:]/fcor,
        vmin=-0.5,vmax=0.5,shading='auto',cmap=plt.cm.get_cmap('BrBG_r',16), 
        transform=ccrs.PlateCarree(),rasterized=True);                            
# add colorbar                                                             
lon = np.arange(-100,-20,20)
lat = np.array([10,20,30,40,50])
ax.set_xticks(lon, crs=ccrs.PlateCarree())
ax.set_yticks(lat, crs=ccrs.PlateCarree())
lon_formatter = LongitudeFormatter()
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)

axpos = ax.get_position()                                                 
cbar_ax = fig.add_axes([axpos.x1+0.01,axpos.y0+0.005,0.03,axpos.height*0.98]);
cbar = fig.colorbar(im1, cax=cbar_ax)                             
cbar.ax.tick_params(labelsize=16)                                          
cbar.set_label(r'$\xi$/f',rotation=0, fontsize=16)                      

plt.savefig('paperfigs/SST_vorticity.png', bbox_inches='tight',format='png',dpi=300)

