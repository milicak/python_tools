import numpy as np
import cartopy.crs as ccrs
import cartopy.feature
import mpl_util
from cpttoseg import cpt2seg


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



dfs = xr.open_dataset('/archive/milicak/MITgcm_c65/Projects/Arctic_4km/ncfiles/2DArcticOcean_1999_06.nc')
df = xr.open_dataset('/archive/milicak/MITgcm_c65/Projects/Arctic_4km/Exp02_0/2DArcticOcean_2009_07.nc')
lon = np.copy(dfs.XC)
lat = np.copy(dfs.YC)
fcor = coriolis(lat)

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
im1 = m.pcolormesh(lon,lat,ma.masked_where(dfs.Depth==0,df.momVort3[-1,:,:]/fcor),
                   cmap=plt.cm.get_cmap('BrBG_r',16),vmin=-.5,vmax=.5,latlon=True, rasterized=True)
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
im2 = m.pcolormesh(lon,lat,ma.masked_where(dfs.Depth==0,df.THETA[-1,:,:]),
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


plt.savefig('paperfigs/ctrl_vorticity_sst_snap_BrBG.png', bbox_inches='tight',format='png',dpi=300)


# fig = plt.figure(figsize=(6,6))
# im1 = m.pcolormesh(lon,lat,ma.masked_where(dfs.Depth==0,dfs.momVort3[6,:,:]/fcor),
#                    cmap=plt.cm.get_cmap('RdBu_r',16),vmin=-.5,vmax=.5,latlon=True, rasterized=True)
# cb = m.colorbar(im1,"right", size="5%", pad="4%")
# cb.ax.tick_params(labelsize=14)
# parallels = np.arange(40.,86,5.)
# # labels = [left,right,top,bottom]
# m.drawparallels(parallels,labels=[True,False,False,False],fontsize=14);
# meridians = np.arange(10.,351.,20.)
# m.drawmeridians(meridians,labels=[False,False,False,True],fontsize=14);
#
# plt.savefig('paperfigs/ctrl_vorticity_snap.png', bbox_inches='tight',format='png',dpi=300)
