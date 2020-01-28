import numpy as np
#%matplotlib inline
#np.shape !!!!!
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import scipy.io
import pylab
import numpy.ma as ma
from disc_cb import discrete_cmap
#import my_nanfilter
#from my_nanfilter import my_nanfilterbox
import nccf
from netCDF4 import Dataset
import sys
from netcdf_functions import nc_read
from netcdf_functions import ncgetdim
plt.ion()

nx=360
ny=200
nz=50

pylab.ion()
def ncread_time(fname,variable,timestr,timeend,x,y):
    # how to use this subroutine is from netcdf_functions import nc_read
    ncfile = Dataset(fname, 'r', format='NETCDF4')
    tmp=np.zeros([y,x])
    for i in range(timestr,timeend):
        #print i
        tmp=tmp+ncfile.variables[variable][i,0:33,:,:].copy();

    tmp=tmp/(timeend-timestr)
    return tmp

#mask_file='/mnt/fimmhome/Analysis/mom/APE/SO/Analysis/grid_spec_v6_regMask.nc';
#mask=nc_read(mask_file,'tmask');
# Atlantic mask==2

#root_folder='/export/grunchfs/unibjerknes/milicak/bckup/mom/'
root_folder='/shared/projects/uniklima/globclim/milicak/mom/'

projects=['om3_core3_ctrl','om3_core3_patchy_full_01','om3_core3_patchy_full_02']
hist_folder = ['history_63-124years'];

#tempwoa = nc_read('/fimm/home/bjerknes/milicak/Analysis/mom/patchy_NA/Analysis/WOA09_ann_theta_cm2m_extrap.nc','POTENTIAL_TEMP');
tempwoa = nc_read('/shared/projects/uniklima/globclim/milicak/Analysis/mom/patchy_NA/Analysis/WOA09_ann_theta_cm2m_extrap.nc','POTENTIAL_TEMP');
tempwoa = np.squeeze(tempwoa[:,0:33,:,:])
dzwoa = nc_read('/shared/projects/uniklima/globclim/milicak/Analysis/mom/patchy_NA/Analysis/WOA09_ann_theta_cm2m_extrap.nc','DEPTH_bnds');
dzwoa = dzwoa[:,1]-dzwoa[:,0];
dzwoa = dzwoa[0:33]
# instead of repmat use np.tile
dzwoa = np.tile(np.copy(dzwoa),(nx,ny,1))
dzwoa = np.swapaxes(dzwoa,2,0)
lon = nc_read('/shared/projects/uniklima/globclim/milicak/Analysis/mom/patchy_NA/Analysis/levitus_ewg_temp_salt_cm2m.nc','x_T');
lat = nc_read('/shared/projects/uniklima/globclim/milicak/Analysis/mom/patchy_NA/Analysis/levitus_ewg_temp_salt_cm2m.nc','y_T');
tempwoa_dep = tempwoa*dzwoa
tempwoa_dep=np.squeeze(np.nansum(tempwoa_dep,0))/793.625;

# compute amoc    
for i in range(0,1):
    filename=root_folder+projects[i]+'/'+hist_folder[0]+'/'+'temp_00630101.ocean_month.nc'            
    filename2=root_folder+projects[i]+'/'+hist_folder[0]+'/'+'dht_00630101.ocean_month.nc'            
    print filename
    temp=ncread_time(filename,'temp',700,744,nx,ny)
    #temp=ncread_time(filename,'temp',348,744,nx,ny)
    dz=ncread_time(filename2,'dht',700,744,nx,ny)
    tempmom_dep=temp*dz
    tempmom_dep=np.squeeze(np.nansum(tempmom_dep,0))/np.squeeze(np.nansum(dz,0))
    if i==0:
        temp_ctrl=np.copy(tempmom_dep)
    if i==2:
        temp2=np.copy(tempmom_dep)


        
    # bias from WOA
    fig = plt.figure()
    m=Basemap(llcrnrlon=-280,llcrnrlat=-88,urcrnrlon=80,urcrnrlat=88,projection='cyl')
    m.drawcoastlines()
    m.fillcontinents()
    m.drawparallels(np.arange(-80,81,20),labels=[1,0,0,0])
    #m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
    m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
    im1 = m.pcolormesh(lon,lat,np.ma.masked_invalid(tempmom_dep-tempwoa_dep),shading='flat',cmap=discrete_cmap(32, 'RdBu_r')
                      ,vmin=-4,vmax=4,latlon=True)
    cb = m.colorbar(im1,"right", size="5%", pad="10%",ticks=[-4, -3, -2, -1, 0, 1, 2, 3, 4]) # pad is the distance between colorbar and figure
    #cb = plt.colorbar(im1,ticks=[-4,-3, -2, -1, 0, 1, 2, 3, 4]) # pad is the distance between colorbar and figure
    #cb = plt.colorbar(im1,fraction=0.06) # pad is the distance between colorbar and figure
    cb.set_label('[' r'$^\circ$' 'C]')
    plt.show() 
    plt.savefig('paperfigs/'+projects[i]+'_temp33_bias.eps', bbox_inches='tight',format='eps', dpi=300)
    #plt.clf()
    #plt.close(fig)

# difference from sst_ctrl
#fig = plt.figure()
#m=Basemap(llcrnrlon=-280,llcrnrlat=-88,urcrnrlon=80,urcrnrlat=88,projection='cyl')
#m.drawcoastlines()
#m.fillcontinents()
#m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
#m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
#im1 = m.pcolormesh(lon,lat,np.ma.masked_invalid(temp2-temp_ctrl),shading='flat',cmap='RdBu_r'
#                    ,vmin=-4,vmax=4,latlon=True)
##cb = m.colorbar(im1,"right", size="5%", pad="10%",ticks=[-4,-3, -2, -1, 0, 1, 2, 3, 4]) # pad is the distance between colorbar and figure
#cb = plt.colorbar(im1,ticks=[-4,-3, -2, -1, 0, 1, 2, 3, 4]) # pad is the distance between colorbar and figure
#cb.set_label('[' r'$^\circ$' 'C]')
#plt.show() 
#plt.savefig('paperfigs/'+projects[i]+'_temp33_diff_patchy2_ctrl.eps', bbox_inches='tight',format='eps', dpi=300)
#plt.clf()
#plt.close(fig)
