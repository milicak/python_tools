import numpy as np
#%matplotlib inline
#np.shape !!!!!
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import scipy.io
import numpy.ma as ma
from disc_cb import discrete_cmap
#import my_nanfilter
#from my_nanfilter import my_nanfilterbox
import nccf
from netCDF4 import Dataset
import sys
from netcdf_functions import nc_read
from netcdf_functions import ncgetdim

nx=360
ny=200
nz=50

def ncread_time_surface(fname,variable,timestr,timeend,x,y):
    # how to use this subroutine is from netcdf_functions import nc_read
    ncfile = Dataset(fname, 'r', format='NETCDF4')
    tmp=np.zeros([y,x])
    for i in range(timestr,timeend):
        #print i
        tmp=tmp+ncfile.variables[variable][i,0,:,:].copy();

    tmp=tmp/(timeend-timestr)
    return tmp

#mask_file='/mnt/fimmhome/Analysis/mom/APE/SO/Analysis/grid_spec_v6_regMask.nc';
#mask=nc_read(mask_file,'tmask');
# Atlantic mask==2

root_folder='/mnt/fimm/mom/'

projects=['om3_core3_ctrl','om3_core3_patchy_full_01','om3_core3_patchy_full_02']
hist_folder = ['history_63-124years'];

legendnames=['Cold-ctrl','Exp1','Exp2','Exp3','Warm-ctrl','Exp4','Exp5']
plotcolors=['cyan','blue','magenta','brown','green','red','black']

saltwoa = nc_read('/mnt/fimmhome/Analysis/mom/patchy_NA/Analysis/WOA09_ann_salinity_cm2m_extrap.nc','S_AN');
ssswoa=np.squeeze(saltwoa[:,0,:,:])
lon = nc_read('/mnt/fimmhome/Analysis/mom/patchy_NA/Analysis/levitus_ewg_temp_salt_cm2m.nc','x_T');
lat = nc_read('/mnt/fimmhome/Analysis/mom/patchy_NA/Analysis/levitus_ewg_temp_salt_cm2m.nc','y_T');

# compute amoc    
for i in range(0,3):
    filename=root_folder+projects[i]+'/'+hist_folder[0]+'/'+'00630101.ocean_month.nc'            
    print filename
    #sss=ncread_time_surface(filename,'salt',740,744,nx,ny)
    sss=ncread_time_surface(filename,'salt',348,744,nx,ny)
    if i==0:
        sss_ctrl=np.copy(sss)
    if i==2:
        sss2=np.copy(sss)
        
        
    #sss=ncread_time_surface(filename,'salt',740,744,nx,ny)
    # bias from WOA
    fig = plt.figure()
    m=Basemap(llcrnrlon=-280,llcrnrlat=-88,urcrnrlon=80,urcrnrlat=88,projection='cyl')
    m.drawcoastlines()
    m.fillcontinents()
    m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
    m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
    im1 = m.pcolormesh(lon,lat,np.ma.masked_invalid(sss-ssswoa),shading='flat',cmap='RdBu_r' #cmap=discrete_cmap(11, 'RdBu_r')
                      ,vmin=-3,vmax=3,latlon=True)
    cb = m.colorbar(im1,"right", size="5%", pad="10%",ticks=[-3, -2, -1, 0, 1, 2, 3]) # pad is the distance between colorbar and figure
    cb.set_label('[psu]')
    plt.show() 
    plt.savefig('paperfigs/'+projects[i]+'_sss_bias.eps', bbox_inches='tight',format='eps', dpi=300)
    plt.clf()
    plt.close(fig)

# difference from sss_ctrl
fig = plt.figure()
m=Basemap(llcrnrlon=-280,llcrnrlat=-88,urcrnrlon=80,urcrnrlat=88,projection='cyl')
m.drawcoastlines()
m.fillcontinents()
m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
im1 = m.pcolormesh(lon,lat,np.ma.masked_invalid(sss2-sss_ctrl),shading='flat',cmap='RdBu_r'
                    ,vmin=-2,vmax=2,latlon=True)
cb = m.colorbar(im1,"right", size="5%", pad="10%",ticks=[-3, -2, -1, 0, 1, 2, 3]) # pad is the distance between colorbar and figure
cb.set_label('[psu]')
plt.show() 
plt.savefig('paperfigs/'+projects[i]+'_sss_diff_patchy2_ctrl.eps', bbox_inches='tight',format='eps', dpi=300)
plt.clf()
plt.close(fig)
