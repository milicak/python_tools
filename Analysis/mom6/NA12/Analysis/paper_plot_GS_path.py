import numpy as np
import glob
from mpl_toolkits.basemap import Basemap

# GS path: observation [283, 312, 34, 41]
gs_clim = np.loadtxt('GS_clim.txt')
lon0 = -75 # lon0 = 360 - 75
lone = -50 # lone = 360 - 50
nlon = lone - lon0 + 1
gs_lon = np.arange(lon0,lone+1)

root_folder = '/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/'
gr = xr.open_dataset(root_folder + 'ocean_geometry.nc')
df = xr.open_dataset('GS_separation_mom6.nc')

fig, ax1 = plt.subplots(figsize=(8,8))
# m = Basemap(projection='merc',llcrnrlat=15,urcrnrlat=55,
#             llcrnrlon=-100,urcrnrlon=-50,lat_ts=20,resolution='i')
m = Basemap(projection='merc',llcrnrlat=25,urcrnrlat=47,
            llcrnrlon=-90,urcrnrlon=-50,lat_ts=20,resolution='i')
m.drawcoastlines();
m.fillcontinents(color='grey');
im1 = m.pcolormesh(gr.geolon,gr.geolat,df.mean_sst,shading='gouraud',
        vmin=-1,vmax=30,cmap='nice_gfdl',latlon=True);
m.scatter(gs_lon, gs_clim, c='g',latlon=True);
m.contour(gr.geolon,gr.geolat,df.mean_200mtemp,[15],colors='k',linewidths=2,latlon=True)
cb = m.colorbar(im1,"right", size="5%", pad="2%")
cb.ax.tick_params(labelsize=14)
m.contour(gr.geolon[700:900,200:400],gr.geolat[700:900,200:400],
        gr.D[700:900,200:400],levels=[100,300,500,700,1000],colors='gray',linewidths=0.5,latlon=True)

parallels = np.arange(25.,50.,5.)
m.drawparallels(parallels,labels=[True,False,False,False],fontsize=14);
meridians = np.arange(-90.,-45.,10.)
m.drawmeridians(meridians,labels=[False,False,False,True],fontsize=14);

plt.savefig('paperfigs/GS_separation.png', bbox_inches='tight',format='png',dpi=300)


