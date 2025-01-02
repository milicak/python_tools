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
df = xr.open_dataset('ocean_daily_v2.nc')

ssh = np.copy(df.zos)

#Creating the map object
timeind = 29
fig = plt.figure()
m1 = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
m1.set_extent([33.5,36.5,35,37], crs=ccrs.PlateCarree())
m1.coastlines(resolution='10m')
m1.add_feature(cfeature.LAND,color="grey",zorder=2)
m1.add_feature(cfeature.LAKES,zorder=2, alpha=0.5)
m1.add_feature(cfeature.RIVERS,zorder=2)
cm = m1.pcolormesh(gr.geolon,gr.geolat,ssh[timeind,:,:],vmin=-0.25, vmax=.25,
        shading='gouraud',cmap='BrBG_r', transform=ccrs.PlateCarree())
plt.colorbar(cm,fraction=0.03)
m1.contour(gr.geolon,gr.geolat,gr.D,levels=[50,100,150,200,250],colors='grey',
           linewidths=0.25,linestyles='dashed',transform=ccrs.PlateCarree())
plt.title('time = 15 minutes')
fname = 'paperfigs/Antakya_tsunami2_ssh_15mins.png'
plt.savefig(fname, bbox_inches='tight',format='png',dpi=300)

timeind = 59
fig = plt.figure()
m1 = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
m1.set_extent([33.5,36.5,35,37], crs=ccrs.PlateCarree())
m1.coastlines(resolution='10m')
m1.add_feature(cfeature.LAND,color="grey",zorder=2)
m1.add_feature(cfeature.LAKES,zorder=2, alpha=0.5)
m1.add_feature(cfeature.RIVERS,zorder=2)
cm = m1.pcolormesh(gr.geolon,gr.geolat,ssh[timeind,:,:],vmin=-0.25, vmax=.1,
        shading='gouraud',cmap='BrBG_r', transform=ccrs.PlateCarree())
plt.colorbar(cm,fraction=0.03)
m1.contour(gr.geolon,gr.geolat,gr.D,levels=[50,100,150,200,250],colors='grey',
           linewidths=0.25,linestyles='dashed',transform=ccrs.PlateCarree())
plt.title('time = 30 minutes')
fname = 'paperfigs/Antakya_tsunami2_ssh_30mins.png'
plt.savefig(fname, bbox_inches='tight',format='png',dpi=300)

timeind = 119
fig = plt.figure()
m1 = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
m1.set_extent([33.5,36.5,35,37], crs=ccrs.PlateCarree())
m1.coastlines(resolution='10m')
m1.add_feature(cfeature.LAND,color="grey",zorder=2)
m1.add_feature(cfeature.LAKES,zorder=2, alpha=0.5)
m1.add_feature(cfeature.RIVERS,zorder=2)
cm = m1.pcolormesh(gr.geolon,gr.geolat,ssh[timeind,:,:],vmin=-0.25, vmax=.25,
        shading='gouraud',cmap='BrBG_r', transform=ccrs.PlateCarree())
plt.colorbar(cm,fraction=0.03)
m1.contour(gr.geolon,gr.geolat,gr.D,levels=[50,100,150,200,250],colors='grey',
           linewidths=0.25,linestyles='dashed',transform=ccrs.PlateCarree())
plt.title('time = 60 minutes')
fname = 'paperfigs/Antakya_tsunami2_ssh_60mins.png'
plt.savefig(fname, bbox_inches='tight',format='png',dpi=300)

timeind = 239
fig = plt.figure()
m1 = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
m1.set_extent([33.5,36.5,35,37], crs=ccrs.PlateCarree())
m1.coastlines(resolution='10m')
m1.add_feature(cfeature.LAND,color="grey",zorder=2)
m1.add_feature(cfeature.LAKES,zorder=2, alpha=0.5)
m1.add_feature(cfeature.RIVERS,zorder=2)
cm = m1.pcolormesh(gr.geolon,gr.geolat,ssh[timeind,:,:],vmin=-.25, vmax=.25,
        shading='gouraud',cmap='BrBG_r', transform=ccrs.PlateCarree())
plt.colorbar(cm,fraction=0.03)
m1.contour(gr.geolon,gr.geolat,gr.D,levels=[50,100,150,200,250],colors='grey',
           linewidths=0.25,linestyles='dashed',transform=ccrs.PlateCarree())
plt.title('time = 120 minutes')
fname = 'paperfigs/Antakya_tsunami2_ssh_120mins.png'
plt.savefig(fname, bbox_inches='tight',format='png',dpi=300)

