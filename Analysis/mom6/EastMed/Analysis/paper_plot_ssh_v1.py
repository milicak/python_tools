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
df = xr.open_dataset('ocean_daily_v1.nc')

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
cm = m1.pcolormesh(gr.geolon,gr.geolat,ssh[timeind,:,:],vmin=-.1, vmax=.1,
        shading='gouraud',cmap='BrBG_r', transform=ccrs.PlateCarree())
plt.colorbar(cm,fraction=0.03)
m1.contour(gr.geolon,gr.geolat,gr.D,levels=[50,100,150,200,250],colors='grey',
           linewidths=0.25,linestyles='dashed',transform=ccrs.PlateCarree())
plt.title('time = 15 minutes')
fname = 'paperfigs/Antakya_tsunami1_ssh_15mins.png'
plt.savefig(fname, bbox_inches='tight',format='png',dpi=300)

timeind = 59
fig = plt.figure()
m1 = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
m1.set_extent([33.5,36.5,35,37], crs=ccrs.PlateCarree())
m1.coastlines(resolution='10m')
m1.add_feature(cfeature.LAND,color="grey",zorder=2)
m1.add_feature(cfeature.LAKES,zorder=2, alpha=0.5)
m1.add_feature(cfeature.RIVERS,zorder=2)
cm = m1.pcolormesh(gr.geolon,gr.geolat,ssh[timeind,:,:],vmin=-.1, vmax=.1,
        shading='gouraud',cmap='BrBG_r', transform=ccrs.PlateCarree())
plt.colorbar(cm,fraction=0.03)
m1.contour(gr.geolon,gr.geolat,gr.D,levels=[50,100,150,200,250],colors='grey',
           linewidths=0.25,linestyles='dashed',transform=ccrs.PlateCarree())
plt.title('time = 30 minutes')
fname = 'paperfigs/Antakya_tsunami1_ssh_30mins.png'
plt.savefig(fname, bbox_inches='tight',format='png',dpi=300)

timeind = 119
fig = plt.figure()
m1 = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
m1.set_extent([33.5,36.5,35,37], crs=ccrs.PlateCarree())
m1.coastlines(resolution='10m')
m1.add_feature(cfeature.LAND,color="grey",zorder=2)
m1.add_feature(cfeature.LAKES,zorder=2, alpha=0.5)
m1.add_feature(cfeature.RIVERS,zorder=2)
cm = m1.pcolormesh(gr.geolon,gr.geolat,ssh[timeind,:,:],vmin=-.05, vmax=.05,
        shading='gouraud',cmap='BrBG_r', transform=ccrs.PlateCarree())
plt.colorbar(cm,fraction=0.03)
m1.contour(gr.geolon,gr.geolat,gr.D,levels=[50,100,150,200,250],colors='grey',
           linewidths=0.25,linestyles='dashed',transform=ccrs.PlateCarree())
plt.title('time = 60 minutes')
fname = 'paperfigs/Antakya_tsunami1_ssh_60mins.png'
plt.savefig(fname, bbox_inches='tight',format='png',dpi=300)



fig = plt.figure()
m1 = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
m1.set_extent([33.5,36.5,35,37], crs=ccrs.PlateCarree())
m1.coastlines(resolution='10m')
m1.add_feature(cfeature.LAND,color="grey",zorder=2)
m1.add_feature(cfeature.LAKES,zorder=2, alpha=0.5)
m1.add_feature(cfeature.RIVERS,zorder=2)
cm = m1.pcolormesh(gr.geolon,gr.geolat,ssh[timeind,:,:],vmin=-.1, vmax=.1,
        shading='gouraud',cmap='BrBG_r', transform=ccrs.PlateCarree())
plt.colorbar(cm,fraction=0.03)
m1.contour(gr.geolon,gr.geolat,gr.D,levels=[50,100,150,200,250],colors='grey',
           linewidths=0.25,linestyles='dashed',transform=ccrs.PlateCarree())
plt.title('time = 60 minutes')
size = 0.2
iax = plt.axes([0.1, 0.6, size, size], projection=ccrs.PlateCarree(), label='2')
iax.coastlines()
iax.add_feature(cfeature.LAND,color="grey",zorder=2)
iax.pcolormesh(gr.geolon,gr.geolat,ssh[timeind,:,:],vmin=-.1, vmax=.1,shading='gouraud',cmap='BrBG_r', transform=ccrs.PlateCarree())
iax.set_extent([34.73,35.23,36.4,36.76], crs=ccrs.PlateCarree())
m1.indicate_inset_zoom(iax,edgecolor="black");
fname = 'paperfigs/Antakya_tsunami1_ssh_60mins_v2.png'
plt.savefig(fname, bbox_inches='tight',format='png',dpi=300)

timeind = 239
fig = plt.figure()
m1 = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
m1.set_extent([33.5,36.5,35,37], crs=ccrs.PlateCarree())
m1.coastlines(resolution='10m')
m1.add_feature(cfeature.LAND,color="grey",zorder=2)
m1.add_feature(cfeature.LAKES,zorder=2, alpha=0.5)
m1.add_feature(cfeature.RIVERS,zorder=2)
cm = m1.pcolormesh(gr.geolon,gr.geolat,ssh[timeind,:,:],vmin=-.04, vmax=.04,
        shading='gouraud',cmap='BrBG_r', transform=ccrs.PlateCarree())
plt.colorbar(cm,fraction=0.03)
m1.contour(gr.geolon,gr.geolat,gr.D,levels=[50,100,150,200,250],colors='grey',
           linewidths=0.25,linestyles='dashed',transform=ccrs.PlateCarree())
plt.title('time = 120 minutes')
fname = 'paperfigs/Antakya_tsunami1_ssh_120mins.png'
plt.savefig(fname, bbox_inches='tight',format='png',dpi=300)

