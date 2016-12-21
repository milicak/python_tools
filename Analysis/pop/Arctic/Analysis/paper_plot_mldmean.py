import numpy as np
#%matplotlib inline
#np.shape !!!!!
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import scipy.io
#import numpy.ma as ma
import my_nanfilter
from my_nanfilter import my_nanfilterbox
from disc_cb import discrete_cmap
from needJet2 import shfn
from netCDF4 import Dataset
import sys
from netcdf_functions import nc_read
from netcdf_functions import ncgetdim

grid_file='/mnt/fimm/Analysis/pop/Arctic/Analysis/areacello_fx_CCSM4_piControl_r0i0p0.nc';
lon=nc_read(grid_file,'lon');
lat=nc_read(grid_file,'lat');

root_folder='/mnt/fimm/Analysis/pop/Arctic/Analysis/matfiles/'

projects=['B1850CN_f19_tn11_kdsens','B1850CN_f19_tn11_kdsens01','B1850CN_f19_tn11_kdsens02','B1850CN_f19_tn11_kdsens03',
          'B1850CN_f19_tn11_kdsens05','B1850CN_f19_tn11_kdsens04','B1850CN_f19_tn11_kdsens06']

#05 ctrl for warm; 04 Barents,Kara + Eurasia; 06 Barents,Kara for warm exps

av1=75; #240; %75;
av2=85; #250; %85;

fyear = '1'; # first year
lyear = '250'; # last year
parallels = np.arange(-80.,90,20.)
meridians = np.arange(0.,360.,20.)
cmap_needjet2=shfn()

for i in xrange(0,7):
    fig = plt.figure()
    filename=root_folder+projects[i]+'_mld_mean_'+fyear+'_'+lyear+'.mat'
    mat = scipy.io.loadmat(filename)
    mld=np.array(mat['MLDMarch'])
    mld=np.squeeze(np.nanmean(np.squeeze(mld[av1-1:av2-1,:,:]),0))
    mld=mld*1e-2 #convert from cm to meter
    m = Basemap(width=8000000,height=8000000,
            resolution='l',projection='stere',\
            lat_ts=40,lat_0=90,lon_0=0.)
    m.drawcoastlines() #m.drawcoastlines(linewidth=1.5)
    m.fillcontinents()
    m.drawparallels(parallels)
    m.drawmeridians(meridians)
    #im1 = m.pcolormesh(lon,lat,np.transpose(np.ma.masked_invalid(ICEareaMarch)),shading='flat',cmap=plt.cm.gist_ncar_r,latlon=True)
    #im1 = m.pcolormesh(lon,lat,np.transpose(np.ma.masked_invalid(ICEareaMarch)),shading='flat',cmap=plt.cm.spectral,latlon=True)
    #im1 = m.pcolormesh(lon,lat,np.transpose(np.ma.masked_invalid(ICEareaMarch)),shading='flat',cmap=plt.cm.Blues,
    #im1 = m.pcolormesh(lon,lat,np.transpose(np.ma.masked_invalid(ICEareaMarch)),shading='flat',cmap=discrete_cmap(9, 'gist_ncar_r')
    im1 = m.pcolormesh(lon,lat,np.transpose(np.ma.masked_invalid(mld)),shading='flat',cmap=cmap_needjet2
                      ,vmin=0,vmax=1400,latlon=True)
    # vmax=100, vmin=20,latlon=True)
    #cmap = plt.get_cmap('Blues', 10)
    cb = m.colorbar(im1,"right", size="5%", pad="2%",ticks=[0, 200, 400, 600, 800, 1000, 1200, 1400]) # pad is the distance between colorbar and figure
    cb.set_label('[m]')
    #plt.clim(20,110)
    plt.show()
    plt.savefig('paperfigs/'+projects[i]+'_mld_March_'+str(av1)+'_'+str(av2)+'_years.eps', format='eps', dpi=300)    
    plt.clf()
    plt.close(fig)
        
