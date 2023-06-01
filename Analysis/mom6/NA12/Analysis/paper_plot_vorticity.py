import numpy as np
import cartopy.crs as ccrs
import cartopy.feature
import mpl_util
from cpttoseg import cpt2seg
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
df = xr.open_dataset('/archive/milicak/MOM6-examples/Projects/Arctic_Copernicus/20090101.ocean_daily_2009_06.nc',decode_times=False)
gr = xr.open_dataset('/archive/milicak/MOM6-examples/Projects/Arctic_Copernicus/ocean_geometry.nc')
lon = np.copy(gr.geolon)
lat = np.copy(gr.geolat)
lonb = np.copy(gr.geolonb)
latb = np.copy(gr.geolatb)
fcor = coriolis(latb)
gr = gr.rename({'lonh':'xh','lath':'yh','lonq':'xq','latq':'yq'})

# gr2 = xr.open_dataset('/archive/milicak/MOM6-examples/Projects/Arctic_Copernicus/ocean_hgrid.nc')
# angle_dx = gr2.angle_dx[1::2,1::2]

# compute vorticity
df = df.merge(gr)

# 2D grid
grid = Grid(df, coords={'X': {'center': 'xh', 'outer': 'xq'},
                        'Y': {'center': 'yh', 'outer': 'yq'},
                         }, periodic=[])
vorticity = ( - grid.diff(df.ssu * df.dxCu, 'Y', boundary='fill')
              + grid.diff(df.ssv * df.dyCv, 'X', boundary='fill') ) / df.Aq

m = Basemap(width=2600000,height=2300000,
            resolution='l',projection='stere',
            lat_ts=50,lat_0=70,lon_0=0.)
lon1,lat1 = m(-46,74)

fig, axes = plt.subplots(figsize=(9,6))
ax1 = plt.subplot2grid(shape=(1,4),loc=(0,0),colspan=2)
ax2 = plt.subplot2grid(shape=(1,4),loc=(0,2),colspan=2)
plt.tight_layout()

m.ax=ax1
m.drawcoastlines()
m.fillcontinents(color='grey');
im1 = m.pcolormesh(lonb,latb,vorticity[-1,:,:]/fcor,
                   cmap=plt.cm.get_cmap('BrBG_r',16),vmin=-.25,vmax=.25,latlon=True, rasterized=True)
parallels = np.arange(40.,86,5.)
# labels = [left,right,top,bottom]
m.drawparallels(parallels,labels=[True,False,False,False],fontsize=14);
meridians = np.arange(10.,351.,20.)
m.drawmeridians(meridians,labels=[False,False,False,True],fontsize=14);

axpos = ax1.get_position()
cbar_ax = fig.add_axes([axpos.x0-0.1,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(im1, cax=cbar_ax, ticklocation='left')
cbar.ax.tick_params(labelsize=14)
cbar.set_label(r'f/$\zeta$',rotation=0,y=1.02,labelpad=-48,fontsize=14)
ax1.text(lon1,lat1,'a)',fontsize=14)

m.ax=ax2
m.drawcoastlines()
m.fillcontinents(color='grey');
im2 = m.pcolormesh(lon,lat,ma.masked_where(gr.D==0,df.tos[-1,:,:]),
                   cmap='nice_gfdl',vmin=-1.5,vmax=14,latlon=True, rasterized=True)
parallels = np.arange(40.,86,5.)
# labels = [left,right,top,bottom]
m.drawparallels(parallels,labels=[False,False,False,False],fontsize=14);
meridians = np.arange(10.,351.,20.)
m.drawmeridians(meridians,labels=[False,False,False,True],fontsize=14);
axpos = ax2.get_position()
cbar_ax = fig.add_axes([axpos.x1+0.03,axpos.y0,0.03,axpos.height])
cbar = fig.colorbar(im2, cax=cbar_ax, ticklocation='right')
cbar.ax.tick_params(labelsize=14)
cbar.set_label('$^\circ$C',rotation=0,y=1.07,labelpad=-40,fontsize=14)
ax2.text(lon1,lat1,'b)',fontsize=14);


plt.savefig('paperfigs/vorticity_sst_snap_BrBG.png', bbox_inches='tight',format='png',dpi=300)


