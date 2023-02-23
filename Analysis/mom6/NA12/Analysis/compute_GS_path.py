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

ls1 = sorted(glob.glob(root_folder + '*ocean_month_z*'))

# Take the last 20 years
ls1 = ls1[24:]
df = xr.open_mfdataset(ls1)

# z-leves coordingate zl=9 is 200 meter
mtemp = df.thetao[:,9,:,:].mean('time')
ssttemp = df.thetao[:,0,:,:].mean('time')

# ls1 = sorted(glob.glob(root_folder + '*ocean_annual.nc'))
# mtemp = df.thetao[:,36,:,:].mean('time')
# ssttemp = df.thetao[:,0,:,:].mean('time')

ds = mtemp.to_dataset(name='mean_200mtemp')
ds['mean_sst'] = ssttemp
ds.to_netcdf('GS_separation_mom6.nc')

# plt.figure()
# plt.pcolormesh(gr.geolon,gr.geolat,ssttemp,shading='gouraud',cmap='nice_gfdl');plt.colorbar()
# # plt.axis([282-360,310-360,32,47])
# plt.axis([-90,310-360,28,47])
# plt.scatter(gs_lon, gs_clim, c='g');
# plt.contour(gr.geolon,gr.geolat,mtemp,[15],colors='r',linewidths=2)
# plt.contour(gr.geolon[700:900,200:400],gr.geolat[700:900,200:400],
#         gr.D[700:900,200:400],levels=[100,300,500,700,1000],colors='gray',linewidths=0.5)
#
#
# m = Basemap(projection='merc',llcrnrlat=15,urcrnrlat=55,
#             llcrnrlon=-100,urcrnrlon=-50,lat_ts=20,resolution='i')
