''' computes sea ice area'''
import my_nanfilter
import numpy as np
import numpy.ma as ma
import scipy.io
import sys
#%matplotlib inline
#np.shape !!!!!
from mpl_toolkits.basemap import Basemap
from cpttoseg import cpt2seg
from netCDF4 import Dataset
from netcdf_functions import nc_read
from netcdf_functions import ncgetdim
import matplotlib.pyplot as plt
import matplotlib as mpllib

# IMPORTANT
plt.ion()

_palette_data = cpt2seg('/fimm/home/bjerknes/milicak/python_tools/Analysis/cpt_files/sst.cpt')
palette = mpllib.colors.LinearSegmentedColormap('palette', _palette_data, 256)


def my_nanfilterbox(y, dn):
    import numpy as np
    # box filter
    filtermi = np.ones(dn)
    # bar filter
    #filtermi = np.barlett(dn)
    ny=np.max(np.size(y))-1
    dummy=np.NaN*np.ones(ny+dn)
    if np.float(dn)/2 == np.round(np.float(dn)/2):
       dn0 = dn/2-1
    else:
       dn0 = (dn-1)/2

    dummy[0+dn0:ny+dn0+1] = y
    yfilter=np.ones(ny+1)
    for n in range(0,ny+1):
      fy = dummy[0+n:dn+n]*filtermi
      II = (~np.isnan(fy))
      wt0 = 1/np.mean(filtermi[II])
      yfilter[n] = wt0*np.mean(fy[II])

    return yfilter


fyear = 33; # first year 1980
lyear = 62; # last year 2009
mw = np.array([31,28,31,30,31,30,31,31,30,31,30,31],dtype=np.float)
mw = mw/sum(mw)
nx = 360
ny = 384
grid_file = '/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid.nc'
area = nc_read(grid_file,'parea')
area = area[:-1,:]
root_folder='/work/milicak/mnt/viljework/archive/'
#root_folder='/work/milicak/mnt/norstore/NS2345K/noresm/cases/'
projects = ['NOIIA_T62_tn11_FAMOS_BG_CTR','NOIIA_T62_tn11_FAMOS_BG_POS','NOIIA_T62_tn11_FAMOS_BG_NEG']

ice_area = {}
ice_volume = {}


for project in projects:
    ice_area[project] = []
    ice_volume[project] = []
    for year in xrange(np.int(fyear),np.int(lyear)+1):
        aice = np.zeros([ny,nx])
        vice = np.zeros([ny,nx])
        for month in xrange(1,13):
            filename = root_folder+project+'/ice/hist/'+project+'.cice.h.'+str(year).zfill(4)+'-'+str(month).zfill(2)+'.nc'
            dnm = np.squeeze(nc_read(filename,'aice'))/100.0 # ratio
            # if there is area criteria
            #  dnm[dnm<0.15] = np.nan
            aice = aice+dnm*mw[month-1]*area
            #  ice thickness
            dnm = np.squeeze(nc_read(filename,'hi'))
            #  if there is height criteria
            #  dnm[dnm<0.15] = np.nan
            vice = vice+dnm*mw[month-1]*area
            # get rid off southern hemisphere
            aice[:200,:] = 0.0
            vice[:200,:] = 0.0


        print year
        ice_area[project] = np.append(ice_area[project],aice.sum())
        ice_volume[project] = np.append(ice_volume[project],vice.sum())






