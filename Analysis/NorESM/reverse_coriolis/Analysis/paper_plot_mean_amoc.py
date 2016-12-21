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
regions=['Global', 'Pacific','Atlantic']
    
for i in xrange(0,3):
    filename=root_folder+regions[i]+projects[0]+'_moc_mean_years_'+fyear+'_'+lyear+'.mat'
    mat = scipy.io.loadmat(filename)
    moc=np.array(mat['MOC_mean'])
    depth=np.array(mat['depth'])
    #lat=np.array(mat['lat'])
    lat=np.linspace(-80,85,166)
    moc=np.copy(np.mean(moc,2))
    fig = plt.figure()
    ax = plt.gca()
    ax.set_axis_bgcolor('grey')
    im1 = plt.pcolor(lat,-depth,np.transpose(np.ma.masked_invalid(moc)),cmap='jet',vmin=-16,vmax=28)
    #cb = m.colorbar(im1,"right", size="5%", pad="10%") #,ticks=[-4, -3, -2, -1, 0, 1, 2, 3, 4]) #dnm pad is the distance between colorbar and figure
    cb = plt.colorbar(im1,pad=0.05) 
    cb.set_label('[Sv]')
    plt.ylabel('Depth [m]')
    plt.xlabel('Lat')
    if(i==0): 
      plt.xlim(-80,80)
    else :
      plt.xlim(-40,80)


    plt.ylim(-7000,0)
    #cb.set_label('[' r'$^\circ$' 'C]')
    plt.savefig('paperfigs/'+projects[0]+'_'+regions[i]+'_moc_mean.eps', bbox_inches='tight',format='eps', dpi=200)
    plt.clf()
    plt.close(fig)


#plt.show()



