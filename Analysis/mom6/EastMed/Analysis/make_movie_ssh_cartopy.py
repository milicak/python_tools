#importing required libraries
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pyplot as plt
import cartopy.mpl.geoaxes
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.animation import FuncAnimation
from matplotlib import pyplot as plt, animation

gr = xr.open_dataset('ocean_geometry.nc')
df = xr.open_dataset('ocean_daily.nc')

ssh = np.zeros((df.time.shape[0],df.yh.shape[0],df.xh.shape[0]))
ssh[:,:,:] = np.copy(df.zos)

img = plt.imread('/ari/users/milicak/MTS_logo1.jpg')

#Creating the map object
fig = plt.figure()

for timeind in range(0,df.time.shape[0]):
    print(timeind)
    m1 = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    m1.set_extent([33.5,36.5,35,37], crs=ccrs.PlateCarree())
    m1.coastlines(resolution='10m')
    #Adding features to the map
    m1.add_feature(cfeature.LAND,color="grey",zorder=2)
    #m1.add_feature(cfeature.COASTLINE)
    m1.add_feature(cfeature.BORDERS, zorder=2, linestyle=':')
    m1.add_feature(cfeature.LAKES,zorder=2, alpha=0.5)
    m1.add_feature(cfeature.RIVERS,zorder=2)
    #m1.stock_img()
    cm = m1.pcolormesh(gr.geolon,gr.geolat,ssh[timeind,:,:],vmin=-.25, vmax=.25,
            shading='gouraud',cmap='BrBG_r', transform=ccrs.PlateCarree())
    plt.colorbar(cm,fraction=0.03)
    m1.contour(gr.geolon,gr.geolat,gr.D,levels=[50,100,150,200,250],colors='grey',
               linewidths=0.25,linestyles='dashed',transform=ccrs.PlateCarree())
    m1.imshow(img, origin='upper', extent=(36, 36.5, 35.0, 35.3),zorder=10,transform=ccrs.PlateCarree())
    title = str(df.time[timeind].coords)
    plt.title('Sea surface height ' + title[-5:] + ' minutes')
    fname = 'gifs/Antakya_tsunami' + str(timeind).zfill(3) + '.png'
    plt.savefig(fname, bbox_inches='tight',format='png',dpi=300)
    plt.clf()




# def animate(i):
#     m1.pcolormesh(gr.geolon,gr.geolat,ssh[i,:,:],vmin=-.25, vmax=.25,shading='gouraud',cmap='BrBG_r', transform=ccrs.PlateCarree())
#     m1.imshow(img, origin='upper', extent=(36, 36.5, 35.0, 35.3),zorder=10,transform=ccrs.PlateCarree())
#
#
# anim = animation.FuncAnimation(fig, animate, interval=140, frames=110) # , repeat = False)
# anim.save('Antakya_Tsunami_v2.gif')
# # BEX4_tEJB3.U745
