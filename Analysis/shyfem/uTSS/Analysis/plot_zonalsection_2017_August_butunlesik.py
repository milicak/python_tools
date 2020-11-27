import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import ESMF
from mpl_toolkits.basemap import Basemap                                            
from scipy.interpolate import interp1d
import geopy.distance

plt.ion()

root_folder  = '/work/mi19918/Projects/'
project_name = 'uTSS'

# expid = 'Exp01.2'
# expid = 'Exp_20160101'
# expid = 'Exp_2016_analysis'
# expid = 'Exp_2016_analysis_keps'
expid = 'Exp_2016_analysis_newTSIC'


inds = np.linspace(585,593,593-585+1,dtype=int)

for x,ind in enumerate(inds):
    sdate = "%4.4d" % (ind)
    print(x)
    fname = root_folder+project_name+'/'+expid+'/OUT_2017_2018_2019/'+'uTSS_lobc_chunk_'+sdate+'.nos.nc'
    # fname = root_folder+project_name+'/'+expid+'/OUT/2018/01/'+'uTSS_lobc_chunk_'+sdate+'.nos.nc'
    if x == 0:
        list=glob.glob(fname)
    else:
        list.append(fname)



df = xr.open_mfdataset(list)
ds = df.salinity.mean('time')
dt = df.salinity.mean('time')

m = Basemap(llcrnrlon=26.0,llcrnrlat=40,urcrnrlon=30.,urcrnrlat=41.5,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='h',projection='merc',\
            lat_0=40.,lon_0=20.,lat_ts=20.)

m.drawcoastlines(linewidth=0.2);
m.fillcontinents(color='grey');
m.drawparallels(np.arange(38,44,0.5),labels=[1,1,0,0]);
m.drawmeridians(np.arange(22,33,1),labels=[0,0,0,1]);

longitude,latitude = m(np.copy(df.longitude[-1,:]),np.copy(df.latitude[-1,:]))

im1 = plt.tripcolor(longitude,latitude,df.element_index[-1,:]-1,ds[:,0]
                ,vmin=18,vmax=25
                ,cmap='needJet2',shading='gouraud')

cb = m.colorbar(im1,"right", size="5%", pad="10%")
cb.set_label('$psu$',rotation=0,y=1.07,labelpad=-45)

# plt.savefig('paperfigs/uTSS_SSS.eps', bbox_inches='tight',format='eps', dpi=300)

plt.savefig('paperfigs/uTSS_SSS_2017_08_butunlesik.png',
            bbox_inches='tight',format='png', dpi=300)
