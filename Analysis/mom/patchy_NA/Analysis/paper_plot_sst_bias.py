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
plt.ion()


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

#root_folder='/hexagon/work/milicak/RUNS/mom/om3_core3_2/'
#root_folder='/mnt/fimm/mom/'
root_folder='/export/grunchfs/unibjerknes/milicak/bckup/mom/'

projects=['om3_core3_ctrl','om3_core3_patchy_full_01','om3_core3_patchy_full_02','om3_core3_patchy_full_03']
hist_folder = ['history_63-124years'];
#projects=['om3_core3_patchy_full_03']
#hist_folder = ['history_63-124years'];
#hist_folder = ['history'];

fname1 = '/export/grunchfs/unibjerknes/milicak/bckup/Analysis/mom/patchy_NA/Analysis/WOA09_ann_theta_cm2m_extrap.nc'
fname2 = "/export/grunchfs/unibjerknes/milicak/bckup/Analysis/mom/patchy_NA/Analysis/levitus_ewg_temp_salt_cm2m.nc"
#fname1 = ['/mnt/fimmhome/Analysis/mom/patchy_NA/Analysis/WOA09_ann_theta_cm2m_extrap.nc']
#fname2 = ['/mnt/fimmhome/Analysis/mom/patchy_NA/Analysis/levitus_ewg_temp_salt_cm2m.nc']
tempwoa = nc_read(fname1,'POTENTIAL_TEMP');
lon = nc_read(fname2,'x_T');
lat = nc_read(fname2,'y_T');
sstwoa=np.squeeze(tempwoa[:,0,:,:])

#for i in range(0,1):
for i in range(0,4):
    #filename=root_folder+projects[i]+'/'+hist_folder[0]+'/'+'00010101.ocean_month.nc'
    filename=root_folder+projects[i]+'/'+hist_folder[0]+'/'+'00630101.ocean_month.nc'
    print filename
    #sst=ncread_time_surface(filename,'temp',700,744,nx,ny)
    sst=ncread_time_surface(filename,'temp',348,744,nx,ny)
    if i==0:
        sst_ctrl=np.copy(sst)
    if i==2:
        sst2=np.copy(sst)
    if i==3:
        sst3=np.copy(sst)


    # bias from WOA
    fig = plt.figure()
    m=Basemap(llcrnrlon=-280,llcrnrlat=-88,urcrnrlon=80,urcrnrlat=88,projection='cyl')
    m.drawcoastlines()
    m.fillcontinents()
    m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
    m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
    im1 = m.pcolormesh(lon,lat,np.ma.masked_invalid(sst-sstwoa),shading='flat',cmap=discrete_cmap(32, 'RdBu_r')
                      ,vmin=-4,vmax=4,latlon=True)
    cb = m.colorbar(im1,"right", size="5%", pad="10%",ticks=[-4, -3, -2, -1, 0, 1, 2, 3, 4]) # pad is the distance between colorbar and figure
    cb.set_label('[' r'$^\circ$' 'C]')
    plt.show()
    plt.savefig('paperfigs/'+projects[i]+'_sst_bias.eps', bbox_inches='tight',format='eps', dpi=300)
    plt.clf()
    plt.close(fig)

#sys.exit()
i=2
# difference from sst_ctrl
fig = plt.figure()
m=Basemap(llcrnrlon=-280,llcrnrlat=-88,urcrnrlon=80,urcrnrlat=88,projection='cyl')
m.drawcoastlines()
m.fillcontinents()
m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
im1 = m.pcolormesh(lon,lat,np.ma.masked_invalid(sst2-sst_ctrl),shading='flat',cmap='RdBu_r'
                    ,vmin=-4,vmax=4,latlon=True)
cb = m.colorbar(im1,"right", size="5%", pad="10%",ticks=[-4,-3, -2, -1, 0, 1, 2, 3, 4]) # pad is the distance between colorbar and figure
cb.set_label('[' r'$^\circ$' 'C]')
plt.show()
plt.savefig('paperfigs/'+projects[i]+'_sst_diff_patchy2_ctrl.eps', bbox_inches='tight',format='eps', dpi=300)
plt.clf()
plt.close(fig)

fig = plt.figure()
i=3
m=Basemap(llcrnrlon=-280,llcrnrlat=-88,urcrnrlon=80,urcrnrlat=88,projection='cyl')
m.drawcoastlines()
m.fillcontinents()
m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
im1 = m.pcolormesh(lon,lat,np.ma.masked_invalid(sst3-sst_ctrl),shading='flat',cmap='RdBu_r'
                    ,vmin=-4,vmax=4,latlon=True)
cb = m.colorbar(im1,"right", size="5%", pad="10%",ticks=[-4,-3, -2, -1, 0, 1, 2, 3, 4]) # pad is the distance between colorbar and figure
cb.set_label('[' r'$^\circ$' 'C]')
plt.show()
plt.savefig('paperfigs/'+projects[i]+'_sst_diff_patchy3_ctrl.eps', bbox_inches='tight',format='eps', dpi=300)
plt.clf()
plt.close(fig)
