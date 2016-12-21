import my_nanfilter
import numpy as np
import numpy.ma as ma
import scipy.io
import sys
#%matplotlib inline
#np.shape !!!!!
from needJet2 import shfn
from mpl_toolkits.basemap import Basemap
from cpttoseg import cpt2seg
from my_nanfilter import my_nanfilterbox
from netCDF4 import Dataset
from netcdf_functions import nc_read
from netcdf_functions import ncgetdim
import matplotlib.pyplot as plt
import matplotlib as mpllib
reload(my_nanfilter)

lon=nc_read('/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid.nc','plon')
lat=nc_read('/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid.nc','plat')
lon=np.copy(lon[:-1,:])
lat=np.copy(lat[:-1,:])

cmap_needjet2=shfn()

_palette_data = cpt2seg('/fimm/home/bjerknes/milicak/python_tools/Analysis/cpt_files/sst.cpt')
palette = mpllib.colors.LinearSegmentedColormap('palette', _palette_data, 256)

def enable_global(tlon,tlat,data):
  """Fix the data in such a way that it can to be plotted on a global projection on its native grid"""
  tlon = np.where(np.greater_equal(tlon,min(tlon[:,0])),tlon-360,tlon)
  tlon=tlon+abs(ma.max(tlon)); tlon=tlon+360
  # stack grids side-by-side (in longitiudinal direction), so
  # any range of longitudes may be plotted on a world map.
  tlon = np.concatenate((tlon,tlon+360),1)
  tlat = np.concatenate((tlat,tlat),1)
  data = ma.concatenate((data,data),1)
  tlon = tlon-360.
  return tlon, tlat, data


fyear = '650'; # first year
lyear = '750'; # last year

root_folder='/fimm/home/bjerknes/milicak/Analysis/NorESM/reverse_coriolis/Analysis/matfiles/'

#projects=['N1850_f19_tn11_reverseCoriolis']
projects=['N1850_f19_tn11_01_default']

legendnames=['Cold-ctrl','Exp1','Exp2','Exp3','Warm-ctrl','Exp4','Exp5']
plotcolors=['cyan','blue','magenta','brown','green','red','black']
    
for i in xrange(0,1):
    filename=root_folder+projects[i]+'_timemean_'+fyear+'_'+lyear+'.mat'
    mat = scipy.io.loadmat(filename)
    temp=np.array(mat['salnlvl'])
    sss=np.copy(np.squeeze(temp[:,:-1,0]))
    [lon,lat,sss]=enable_global(lon,lat,np.transpose(sss))
    fig = plt.figure()
    m=Basemap(llcrnrlon=-180,llcrnrlat=-80,urcrnrlon=180,urcrnrlat=90,projection='cyl')
    m.drawcoastlines()
    m.fillcontinents()
    m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
    m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
    im1 = m.pcolormesh(np.transpose(lon-110),np.transpose(lat),np.transpose(np.ma.masked_invalid(sss))
                      ,shading='flat',cmap=cmap_needjet2,vmin=31,vmax=38,latlon=True)
    cb = m.colorbar(im1,"right", size="5%", pad="15%") #,ticks=[-4, -3, -2, -1, 0, 1, 2, 3, 4]) # pad is the distance between colorbar and figure
    cb.set_label('[' r'$^\circ$' 'C]')


#plt.show()
plt.savefig('paperfigs/'+projects[i]+'_sss.eps', bbox_inches='tight',format='eps', dpi=200)
plt.clf()
plt.close(fig)



