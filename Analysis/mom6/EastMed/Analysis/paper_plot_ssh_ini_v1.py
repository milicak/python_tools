import xarray as xr
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt
import cartopy.mpl.geoaxes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# gr = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/EastMed/OUT/ocean_geometry.nc')
# df = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/EastMed/OUT/ocean_daily_v1.nc')

# to plot in altay
gr = xr.open_dataset('ocean_geometry.nc')
df = xr.open_dataset('INPUT/ssh_IC_EastMed.nc')

ssh = np.copy(df.ssh)

#Creating the map object
fig  = plt.figure()
m1 = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
m1.set_extent([33.5,36.5,35,37], crs=ccrs.PlateCarree())
m1.coastlines(resolution='10m')
#Adding features to the map
m1.add_feature(cfeature.LAND,color="grey",zorder=2)
#m1.add_feature(cfeature.COASTLINE)
# m1.add_feature(cfeature.BORDERS, zorder=2, linestyle=':')
m1.add_feature(cfeature.LAKES,zorder=2, alpha=0.5)
m1.add_feature(cfeature.RIVERS,zorder=2)
#m1.stock_img()
timeind = 0
cm = m1.pcolormesh(gr.geolon,gr.geolat,ssh[timeind,:,:],vmin=-.25, vmax=.25,
        shading='gouraud',cmap='BrBG_r', transform=ccrs.PlateCarree())
plt.colorbar(cm,fraction=0.03)
m1.contour(gr.geolon,gr.geolat,gr.D,levels=[50,100,150,200,250],colors='grey',
           linewidths=0.25,linestyles='dashed',transform=ccrs.PlateCarree())

size = 0.2
iax = plt.axes([0.1, 0.6, size, size], projection=ccrs.PlateCarree(), label='2')
iax.coastlines()
iax.pcolormesh(gr.geolon,gr.geolat,ssh[timeind,:,:],vmin=-.25, vmax=.25,shading='gouraud',cmap='BrBG_r', transform=ccrs.PlateCarree())
iax.set_extent([35.72,35.85,36.2,36.32], crs=ccrs.PlateCarree())
m1.indicate_inset_zoom(iax,edgecolor="black");

fname = 'paperfigs/Antakya_tsunami1_initialssh.png'
plt.savefig(fname, bbox_inches='tight',format='png',dpi=300)

