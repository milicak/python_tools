import numpy as np
#%matplotlib inline
#np.shape !!!!!
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib
import scipy.io
#import numpy.ma as ma
import my_nanfilter
from my_nanfilter import my_nanfilterbox
from disc_cb import discrete_cmap
from needJet2 import shfn
import sys
from cpttoseg import cpt2seg
import cdms2
import mpl_util

from netcdf_functions2 import nc_read

#from netCDF4 import Dataset
#from netcdf_functions import ncgetdim

lon = nc_read('/home/mil021/roms_toolbox/matlab/bath/etopo5.nc','topo_lon')
lat = nc_read('/home/mil021/roms_toolbox/matlab/bath/etopo5.nc','topo_lat')
topo = nc_read('/home/mil021/roms_toolbox/matlab/bath/etopo5.nc','topo')
[lons,lats]=np.meshgrid(lon,lat)

parallels = np.arange(-80.,90,20.)
meridians = np.arange(0.,360.,20.)
cmap_needjet2=shfn()

_palette_data = cpt2seg('/home/mil021/python_tools/Analysis/cpt_files/ibcao.cpt')
palette = matplotlib.colors.LinearSegmentedColormap('palette', _palette_data, 24)

#_palette_data = cpt2seg('/home/mil021/python_tools/Analysis/cpt_files/arctic.cpt')
#palette = matplotlib.colors.LinearSegmentedColormap('palette', _palette_data, 110)

conts=[-5000,-4000,-3000,-2500,-2000,-1500,-1000,-500,-200,-100,-50,-25,-10,0,50,100,200,300,400,500,600,700,800,1000]

fig = plt.figure()
m = Basemap(width=5000000,height=5000000,
         resolution='l',projection='stere',
	 lat_ts=40,lat_0=90,lon_0=0.)
m.drawparallels(parallels)
m.drawmeridians(meridians)

#m.drawcoastlines() #m.drawcoastlines(linewidth=1.5)
    #m.fillcontinents()
#im1 = m.pcolormesh(lon,lat,np.transpose(np.ma.masked_invalid(topo)),shading='flat',cmap=palette,
#im1 = m.pcolormesh(lons,lats,topo,shading='flat',cmap=palette,
#          vmax=1400, vmin=-5000,latlon=True)
im1 = m.contourf(lons,lats,topo,conts,cmap=mpl_util.LevelColormap(conts,cmap=palette),
          vmax=1000, vmin=-5000,latlon=True)
cb = m.colorbar(im1,"right", size="5%", pad="2%")

#cmap = plt.get_cmap('Blues', 10)
    #cb.set_label('[' r'$^\circ$' 'C]')
    #plt.clim(20,110)
plt.savefig('paperfigs/Arctic_topo.eps', format='eps', dpi=300)   
plt.show()
 
    #plt.clf()
    #plt.close(fig)
    

