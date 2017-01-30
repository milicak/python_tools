import my_nanfilter
import numpy as np
import numpy.ma as ma
import scipy.io
import sys
#%matplotlib inline
#np.shape !!!!!
from mpl_toolkits.basemap import Basemap
from cpttoseg import cpt2seg
from my_nanfilter import my_nanfilterbox
from netCDF4 import Dataset
from netcdf_functions import nc_read
from netcdf_functions import ncgetdim
import matplotlib.pyplot as plt
import matplotlib as mpllib
reload(my_nanfilter)
sys.path.insert(0,'/fimm/home/bjerknes/milicak/python_tools/Analysis/NorESM/general/Analysis')
from general_diagnostics import timemean
from general_diagnostics import get_sdate_ini

plt.ion()

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


lon=nc_read('/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid.nc','plon')
lat=nc_read('/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid.nc','plat')
lon=np.copy(lon[:-1,:])
lat=np.copy(lat[:-1,:])

_palette_data = cpt2seg('/fimm/home/bjerknes/milicak/python_tools/Analysis/cpt_files/sst.cpt')
palette = mpllib.colors.LinearSegmentedColormap('palette', _palette_data, 256)

fyear = 650; # first year
lyear = 750; # last year

Nx = 360
Ny = 385

root_folder = '/work/milicak/mnt/norstore/NS2345K/noresm/cases/'
#projects = ['N1850_f19_tn11_01_default']
projects = ['N1850_f19_tn11_01_default','N1850_f19_tn11_reverseCoriolis']
mdl = 'micom'
cmpnt = 'ocn'

for n, project in enumerate(projects):
    mld = np.zeros((Ny,Nx))
    maxmld = np.zeros((Ny,Nx))
    if n==0:
        ext = 'hm'
        m2y = 1
    else:
        ext = 'hy'
        m2y = 0


    prefix,sdate = get_sdate_ini(root_folder, cmpnt, mdl, ext, m2y=m2y,
                                 expid=project)
    mld = timemean(prefix, mld, 'mld', fyear, lyear, m2y=m2y,
                   expid=project)
    maxmld = timemean(prefix, maxmld, 'maxmld', fyear, lyear, m2y=m2y,
                   expid=project)


   # sys.exit()
    #[lon,lat,mld]=enable_global(lon,lat,np.transpose(mld))
    #[lon,lat,mld]=enable_global(lon,lat,mld[:-1,:])
    if n==0:
        [lon,lat,maxmld]=enable_global(lon,lat,maxmld[:-1,:])
    else:
        [lon1,lat1,maxmld]=enable_global(lon,lat,maxmld[:-1,:])


    fig = plt.figure()
    m=Basemap(llcrnrlon=-180,llcrnrlat=-80,urcrnrlon=180,urcrnrlat=90,projection='cyl')
    m.drawcoastlines()
    m.fillcontinents()
    m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
    m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
    im1 = m.pcolormesh(np.transpose(lon-110),np.transpose(lat),
                       np.transpose(np.ma.masked_invalid(maxmld))
                      ,shading='flat',cmap=palette,vmin=0,vmax=750,latlon=True)
    cb = m.colorbar(im1,"right", size="5%", pad="15%") #,ticks=[-4, -3, -2, -1, 0, 1, 2, 3, 4]) # pad is the distance between colorbar and figure
    cb.set_label('[m]')
    #plt.show()
    plt.savefig('paperfigs/'+project+'_mld.eps', bbox_inches='tight',format='eps', dpi=200)
    plt.clf()
    plt.close(fig)



