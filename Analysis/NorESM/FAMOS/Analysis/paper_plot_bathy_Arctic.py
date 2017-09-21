''' computes sea ice area'''
#import my_nanfilter
import numpy as np
import numpy.ma as ma
from mpl_toolkits.basemap import Basemap
import scipy.io
import sys
#%matplotlib inline
#np.shape !!!!!
#from mpl_toolkits.basemap import Basemap
from cpttoseg import cpt2seg
from disc_cb import discrete_cmap
from netCDF4 import Dataset
from netcdf_functions import nc_read
from netcdf_functions import ncgetdim
import matplotlib.pyplot as plt
import matplotlib as mpllib
import mpl_util

# IMPORTANT
plt.ion()

_palette_data = cpt2seg('/fimm/home/bjerknes/milicak/python_tools/Analysis/cpt_files/sst.cpt')
palette = mpllib.colors.LinearSegmentedColormap('palette', _palette_data, 256)

_palette_data = cpt2seg('coldblue.cpt')
palette = mpllib.colors.LinearSegmentedColormap('palette', _palette_data, 11)

lon = nc_read('/work/milicak/etopo5.nc','topo_lon')
lat = nc_read('/work/milicak/etopo5.nc','topo_lat')
topo = nc_read('/work/milicak/etopo5.nc','topo')
lon1 = nc_read('/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid.nc','plon');
lat1 = nc_read('/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid.nc','plat');
lon1 = np.transpose(lon1)
lat1 = np.transpose(lat1)
[lons,lats]=np.meshgrid(lon,lat)
parallels = [70,80,85]
meridians = np.arange(0.,360.,45.)

topo[topo>0] = 0

#conts = [-4000,-3000,-2500,-2000,-1500,-1000,-500,-250,-50,0]
conts = [-5000,-4000,-3000,-2500,-2000,-1500,-1000,-500,-250,-50,0]
fig = plt.figure()
m = Basemap(width=6000000,height=6000000,
         resolution='l',projection='stere',
	 lat_ts=45,lat_0=90,lon_0=0.,round='true')
m.drawparallels(parallels)
m.drawmeridians(meridians)
m.drawcoastlines()
m.fillcontinents(color='1')
m.drawrivers()

#im1 = m.contourf(lons,lats,topo,conts,cmap=mpl_util.LevelColormap(conts,cmap=palette),
#im1 = m.contourf(lons,lats,topo,conts,cmap=discrete_cmap(32, 'RdBu_r'),
#im1 = m.contourf(lons,lats,topo,conts,
im1 = m.contourf(lons,lats,topo,conts,cmap=mpl_util.LevelColormap(conts,cmap=palette),
#im1 = m.contourf(lons,lats,topo,conts,cmap=discrete_cmap(24, 'Blues_r'),
          vmax=0, vmin=-5000,latlon=True)
cb = m.colorbar(im1,"right", size="5%", pad="2%")

# barents_opening
bsolon=[lon1[113-1,352-1],lon1[114-1,351-1],lon1[115-1,350-1],lon1[115-1,349-1],
        lon1[116-1,348-1],lon1[117-1,347-1],
        lon1[117-1,346-1],lon1[118-1,345-1],lon1[119-1,344-1],lon1[119-1,343-1],
        lon1[120-1,342-1],lon1[121-1,341-1],
        lon1[121-1,340-1],lon1[122-1,339-1],lon1[122,338]]
bsolat=[lat1[113-1,352-1],lat1[114-1,351-1],lat1[115-1,350-1],lat1[115-1,349-1],
        lat1[116-1,348-1],lat1[117-1,347-1],
        lat1[117-1,346-1],lat1[118-1,345-1],lat1[119-1,344-1],lat1[119-1,343-1],
        lat1[120-1,342-1],lat1[121-1,341-1],
        lat1[121-1,340-1],lat1[122-1,339-1],lat1[122,338]]
xpt,ypt = m(bsolon,bsolat)
#m.plot(xpt,ypt,'b')  # plot a blue dot there

#grid_file = '/work/milicak/ETOPO2v2g_f4.nc'
#lon = nc_read(grid_file,'x')
#lat = nc_read(grid_file,'y')
#topo = nc_read(grid_file,'z')


