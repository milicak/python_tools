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

def plt_arrays(x, y, title="", label="", color="red", linestyle="solid", linewidth=2):
    axes = fig.add_subplot(111)
    axes.plot(x,y, color=color, linestyle=linestyle, linewidth=linewidth,label=label)
    axes.set_title(title)
    #axes.set_label(label)
    #plt.legend(label)
    axes.grid()
    legend = axes.legend(loc='upper right')


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


fyear = '1'; # first year
lyear = '750'; # last year

time = np.linspace(1,750,750)

root_folder='/fimm/home/bjerknes/milicak/Analysis/NorESM/reverse_coriolis/Analysis/matfiles/'

projects=['N1850_f19_tn11_01_default']
#projects=['N1850_f19_tn11_reverseCoriolis']

legendnames=['Global','Pacific','Atlantic']
#plotcolors=['cyan','blue','magenta','brown','green','red','black']
plotcolors=['black','blue','red','brown']
regions=['Global', 'Pacific','Atlantic']
    
fig = plt.figure()
for i in xrange(0,3):
    filename=root_folder+regions[i]+projects[0]+'_mmflxd_years_'+fyear+'_'+lyear+'.mat'
    mat = scipy.io.loadmat(filename)
    moc=np.copy(np.array(mat['amoc_26_5']))
    plt_arrays(time,np.transpose(moc), title="",label=legendnames[i],color=plotcolors[i])


    
    #ax = plt.gca()
    #ax.set_axis_bgcolor('grey')
plt.ylim(0,40)
plt.ylabel('MOC [Sv]')
plt.xlabel('Time [years]')
#plt.show()
plt.savefig('paperfigs/'+projects[0]+'_moc_time.eps', bbox_inches='tight',format='eps', dpi=200)
#plt.clf()
#plt.close(fig)



