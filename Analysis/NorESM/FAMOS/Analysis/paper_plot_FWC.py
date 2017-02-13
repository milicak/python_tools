''' computes sea ice area'''
import numpy as np
import numpy.ma as ma
import scipy.io
import sys
#%matplotlib inline
#np.shape !!!!!
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
from cpttoseg import cpt2seg
from netcdf_functions import nc_read
#import my_nanfilter
#from netcdf_functions import ncgetdim
import matplotlib.pyplot as plt
import matplotlib as mpllib
import pickle

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


Sref = 34.8 # Sref psu
fyear = 33; # 33 first year 1980
lyear = 62; # 62 last year 2009
mw = np.array([31,28,31,30,31,30,31,31,30,31,30,31],dtype=np.float)
mw = mw/sum(mw)
nx = 360
ny = 385
grid_file = '/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid.nc'
area = nc_read(grid_file,'parea')
#area = area[:-1,:]
#mask = nc_read('NorESM_tnx1v2_arctic_mask.nc','mask')
maskfilename = 'BG_GS_masks.mat'
mat = scipy.io.loadmat(maskfilename)
maskGS =np.transpose(np.array(mat['GS_mask']))
maskBG =np.transpose(np.array(mat['BG_mask']))
root_folder='/work/milicak/mnt/viljework/archive/'
#root_folder='/work/milicak/mnt/norstore/NS2345K/noresm/cases/'
#projects = ['NOIIA_T62_tn11_FAMOS_BG_NEG']
projects = ['NOIIA_T62_tn11_FAMOS_BG_CTR','NOIIA_T62_tn11_FAMOS_BG_POS','NOIIA_T62_tn11_FAMOS_BG_NEG',
           'NOIIA_T62_tn11_FAMOS_BG_POS_90s', 'NOIIA_T62_tn11_FAMOS_BG_NEG_90s']

FWC_BG = {}
FWC_GS = {}

for project in projects:
    FWC_BG[project] = []
    FWC_GS[project] = []
    for year in xrange(np.int(fyear),np.int(lyear)+1):
        fwcyBG = np.zeros([ny,nx])
        fwcyGS = np.zeros([ny,nx])
        for month in xrange(1,13):
            if project == 'NOIIA_T62_tn11_FAMOS_BG_NEG' and year == 42 and month == 8:
                filename = root_folder+project+'/ocn/hist/'+project+'.micom.hm.'+str(year).zfill(4)+'-'+str(month-1).zfill(2)+'.nc'
            else:
                filename = root_folder+project+'/ocn/hist/'+project+'.micom.hm.'+str(year).zfill(4)+'-'+str(month).zfill(2)+'.nc'


            dnm1 = np.squeeze(nc_read(filename,'saln'))
            dnm2 = np.squeeze(nc_read(filename,'dz'))
            s1 = np.copy(dnm1)
            s1[s1<0.0] = np.nan
            smask = s1-Sref
            smask[smask >= 0.0] = 0.0
            smask[smask < 0.0] = 1.0
            smask[np.isnan(smask)] = 0.0
            fwcytmp = (Sref-s1)*dnm2*smask/Sref
            fwcyBG = fwcyBG+np.sum(fwcytmp, axis=0)*maskBG*area*mw[month-1]
            fwcyGS = fwcyGS+np.sum(fwcytmp, axis=0)*maskGS*area*mw[month-1]


        print year
        FWC_BG[project] = np.append(FWC_BG[project],fwcyBG.sum())
        FWC_GS[project] = np.append(FWC_GS[project],fwcyGS.sum())
        #FWC[project] = np.append(FWC[project],fwcy)




pickle.dump( FWC_BG, open( "FWC_BG.p", "wb" ) )

