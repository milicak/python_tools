''' plots a map file '''
import numpy as np
import numpy.ma as ma
import scipy.io
import sys
from mpl_toolkits.basemap import Basemap
from cpttoseg import cpt2seg
from netCDF4 import Dataset
from netcdf_functions import nc_read
from netcdf_functions import ncgetdim
import matplotlib.pyplot as plt
import matplotlib as mpllib

plt.ion()

_palette_data = cpt2seg('/fimm/home/bjerknes/milicak/python_tools/Analysis/cpt_files/sst.cpt')
palette = mpllib.colors.LinearSegmentedColormap('palette', _palette_data, 256)

fig = plt.figure(figsize=(20, 8))
m = Basemap(llcrnrlon=-100,llcrnrlat=-80,urcrnrlon=260,urcrnrlat=90,projection='cyl')
m.drawcoastlines()
m.fillcontinents()
m.drawparallels(np.arange(-80,81,30),labels=[1,1,0,0])
m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
plt.savefig('paperfigs/amoc_sketch.eps', bbox_inches='tight',format='eps',
            dpi=300)



