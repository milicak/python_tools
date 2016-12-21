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

# IMPORTANT 
plt.ion()

_palette_data = cpt2seg('/fimm/home/bjerknes/milicak/python_tools/Analysis/cpt_files/sst.cpt')
palette = mpllib.colors.LinearSegmentedColormap('palette', _palette_data, 256)

def ncread_time_surface(fname,variable,timestr,timeend,x,y):
    # how to use this subroutine is from netcdf_functions import nc_read
    ncfile = Dataset(fname, 'r', format='NETCDF4')
    tmp=np.zeros([y,x])
    for i in range(timestr,timeend):
        #print i
        tmp=tmp+ncfile.variables[variable][i,:,:].copy();

    tmp=tmp/(timeend-timestr)
    return tmp


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
nx = 144;
ny = 96;

root_folder='/work/milicak/mnt/norstore/NS2345K/noresm/cases/'

projects=['N1850_f19_tn11_reverseCoriolis']

prect=np.zeros([ny,nx])
    
for year in xrange(np.int(fyear),np.int(lyear)+1):
    #for month in xrange(1,13):
    #filename = root_folder+projects[0]+'atm/hist/'+projects[0]+'cam.h0.'+str(year).zfill(4)+'-'+str(month).zfill(2)+'.nc'
    filename = root_folder+projects[0]+'/atm/hist/'+projects[0]+'.cam.h1.'+str(year).zfill(4)+'-01-01-00000.nc'
    
    precc = ncread_time_surface(filename,'PRECC',0,365,nx,ny)
    precl = ncread_time_surface(filename,'PRECL',0,365,nx,ny)
    #total precipitation
    prect = prect + (precc+precc)
    print year

prect = prect/(np.float(lyear)-np.float(fyear)+1)
# convert m/s to m/year
prect = prect * (365*86400)

lon = nc_read(filename,'lon')
lat = nc_read(filename,'lat')

#sys.exit()

fig = plt.figure()
#ax = plt.gca()
#ax.set_axis_bgcolor('grey')
m=Basemap(llcrnrlon=-180,llcrnrlat=-88,urcrnrlon=180,urcrnrlat=88,projection='cyl')
m.drawcoastlines()
m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
im1 = m.pcolormesh(lon,lat,np.ma.masked_invalid(prect),shading='flat',cmap='RdBu_r'
                      ,vmin=-0.2,vmax=7,latlon=True)
cb = m.colorbar(im1,"right", size="5%", pad="10%") #,ticks=[-4, -3, -2, -1, 0, 1, 2, 3, 4]) #dnm pad is the distance between colorbar and figure
cb.set_label('[m/year]')
#plt.ylabel('Depth [m]')
#plt.xlabel('Lat')
#plt.ylim(-7000,0)
#cb.set_label('[' r'$^\circ$' 'C]')
plt.savefig('paperfigs/'+projects[0]+'_precip_mean.eps', bbox_inches='tight',format='eps', dpi=200)
plt.clf()
plt.close(fig)


#plt.show()



