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

plt.ion()
grid_file='/bcmhsm/milicak/RUNS/noresm/CORE2/Arctic/DATA/ncar-pop/areacello_fx_CCSM4_piControl_r0i0p0.nc';
lon=nc_read(grid_file,'lon');
lat=nc_read(grid_file,'lat');

#root_folder='/mnt/fimm/Analysis/pop/Arctic/Analysis/matfiles/'
root_folder='/export/grunchfs/unibjerknes/milicak/bckup/Analysis/pop/Arctic/Analysis/matfiles/'

projects=['B1850CN_f19_tn11_kdsens','B1850CN_f19_tn11_kdsens01','B1850CN_f19_tn11_kdsens02',
          'B1850CN_f19_tn11_kdsens03','B1850CN_f19_tn11_kdsens05','B1850CN_f19_tn11_kdsens04',
          'B1850CN_f19_tn11_kdsens06','B1850CN_f19_tn11_kdsens07']

#05 ctrl for warm; 04 Barents,Kara + Eurasia; 06 Barents,Kara for warm exps

av1=240; #240; %75;
av2=250; #250; %85;

fyear = '1'; # first year
lyear = '250'; # last year
parallels = np.arange(-80.,90,20.)
meridians = np.arange(0.,360.,20.)
cmap_needjet2=shfn()

for i in xrange(0,8):
    fig = plt.figure()
    filename=root_folder+projects[i]+'_airslp_mean_'+fyear+'_'+lyear+'.mat'
    mat = scipy.io.loadmat(filename)
    slp=np.array(mat['AIRslpmeanwinter'])
    slp=np.squeeze(np.nanmean(np.squeeze(slp[av1-1:av2-1,:,:]),0))
    slp=slp*1e-2 #convert from Pa to hPa
    lon=np.array(mat['lon'])
    lat=np.array(mat['lat'])    
    [lon,lat]=np.meshgrid(lon,lat)
    if i==0:
      slp_ctrl1=np.copy(slp) 
      tmp1=np.copy(slp_ctrl1)
      tmp1[np.isnan(tmp1)]=0
    elif i==4:
      slp_ctrl2=np.copy(slp) 
      tmp2=np.copy(slp_ctrl2)
      tmp2[np.isnan(tmp2)]=0 

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
    im1 = m.pcolormesh(lon,lat,np.transpose(np.ma.masked_invalid(slp))
                        , shading='flat',cmap=cmap_needjet2
                        , vmin=980,vmax=1030,latlon=True)
    # vmax=100, vmin=20,latlon=True)
    #cmap = plt.get_cmap('Blues', 10)
    cb = m.colorbar(im1,"right", size="5%", pad="2%",
                    ticks=[980, 985, 990, 995, 1000, 1005, 1010, 1015, 1020, 1025, 1030]) # pad is the distance between colorbar and figure
    cb.set_label('[hPa]')
    #plt.clim(20,110)
    plt.show()
    plt.savefig('paperfigs/'+projects[i]+'_slp_winter_'
                +str(av1)+'_'+str(av2)+'_years.eps', 
                bbox_inches='tight', format='eps', dpi=300)    
    plt.clf()
    plt.close(fig)
        
for i in xrange(0,8):
     if(i!=0 and i<4):                  
         fig = plt.figure()
         filename=root_folder+projects[i]+'_airslp_mean_'+fyear+'_'+lyear+'.mat'
         print 'mehmet', filename
         mat = scipy.io.loadmat(filename)
         slp_winter=np.array(mat['AIRslpmeanwinter'])
         slp_winter=np.squeeze(np.nanmean(np.squeeze(slp_winter[av1-1:av2-1,:,:]),0))
         m = Basemap(width=8000000,height=8000000,
             resolution='l',projection='stere',\
             lat_ts=40,lat_0=90,lon_0=0.)
         m.drawcoastlines() #m.drawcoastlines(linewidth=1.5)
         m.fillcontinents()
         m.drawparallels(parallels)
         m.drawmeridians(meridians)
         slp_winter=slp_winter*1e-2 #convert from Pa to hPa
         tmp=np.copy(slp_winter)
         tmp[np.isnan(tmp)]=0
         dnm=tmp-tmp1
         dnm[dnm==0]=np.nan
        # im1 = m.pcolormesh(lon,lat,np.transpose(np.ma.masked_invalid(dnm)),shading='flat',cmap=cmap_needjet2,latlon=True)
         im1 = m.pcolormesh(lon,lat,np.transpose(np.ma.masked_invalid(dnm)),
                            shading='flat',cmap=discrete_cmap(19, 'RdBu_r'),
                            vmin=-4,vmax=4,latlon=True)
         cb = m.colorbar(im1,"right", size="5%", pad="2%")
         cb.set_label('[hPa]')
         plt.show()
         plt.savefig('paperfigs/'+projects[i]+'_slp_winter_diff_'
                     +str(av1)+'_'+str(av2)+'_years.eps', 
                     bbox_inches='tight', format='eps', dpi=300)    
         plt.clf()
         plt.close(fig)
     elif(i>=5):
        fig = plt.figure()
        filename=root_folder+projects[i]+'_airslp_mean_'+fyear+'_'+lyear+'.mat'
        print 'mehmet2', filename
        mat = scipy.io.loadmat(filename)
        slp_winter=np.array(mat['AIRslpmeanwinter'])
        slp_winter=np.squeeze(np.nanmean(np.squeeze(slp_winter[av1-1:av2-1,:,:]),0))
        slp_winter=slp_winter*1e-2 #convert from Pa to hPa
        m = Basemap(width=8000000,height=8000000,
        resolution='l',projection='stere',\
        lat_ts=40,lat_0=90,lon_0=0.)
        m.drawcoastlines() #m.drawcoastlines(linewidth=1.5)
        m.fillcontinents()
        m.drawparallels(parallels)
        m.drawmeridians(meridians)
        tmp=np.copy(slp_winter)
        tmp[np.isnan(tmp)]=0
        dnm=tmp-tmp2
        dnm[dnm==0]=np.nan 
        im1 = m.pcolormesh(lon,lat,np.transpose(np.ma.masked_invalid(dnm)),
                           shading='flat',cmap=discrete_cmap(19, 'RdBu_r'),
                           vmin=-4,vmax=4,latlon=True)
        cb = m.colorbar(im1,"right", size="5%", pad="2%")
        cb.set_label('[hPa]')
        plt.show()
        plt.savefig('paperfigs/'+projects[i]+'_slp_winter_diff_'
                    +str(av1)+'_'+str(av2)+'_years.eps', 
                    bbox_inches='tight', format='eps', dpi=300)    
        plt.clf()
        plt.close(fig)
        

