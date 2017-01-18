import numpy as np
#%matplotlib inline
#np.shape !!!!!
import scipy.io as io
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import scipy.io
import numpy.ma as ma
from disc_cb import discrete_cmap
from needJet2 import shfn
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
cmap_needjet2=shfn()

def ncread_time_surface(fname,variable,timestr,timeend,x,y):
    # how to use this subroutine is from netcdf_functions import nc_read
    ncfile = Dataset(fname, 'r', format='NETCDF4')
    tmp=np.zeros([y,x])
    for i in range(timestr,timeend):
        print i
        tmp=tmp+ncfile.variables[variable][i,:,:].copy();

    tmp=tmp/(timeend-timestr)
    return tmp


def ncread_time_surface_ksi(fname,variable,timestr,timeend,x,y,ssh):
    # how to use this subroutine is from netcdf_functions import nc_read
    ncfile = Dataset(fname, 'r', format='NETCDF4')
    tmp=np.zeros([y,x])
    for i in range(timestr,timeend):
        print i
        tmp=tmp+ncfile.variables[variable][i,:,:].copy()-ssh;

    tmp=tmp/(timeend-timestr)
    return tmp


#mask_file='/mnt/fimmhome/Analysis/mom/APE/SO/Analysis/grid_spec_v6_regMask.nc';
#mask=nc_read(mask_file,'tmask');
# Atlantic mask==2

root_folder='/export/grunchfs/unibjerknes/milicak/bckup/mom/'

projects=['om3_core3_ctrl','om3_core3_patchy_full_01','om3_core3_patchy_full_02']
hist_folder = ['history_63-124years'];

lon = nc_read('/fimm/home/bjerknes/milicak/Analysis/mom/patchy_NA/Analysis/levitus_ewg_temp_salt_cm2m.nc','x_T');
lat = nc_read('/fimm/home/bjerknes/milicak/Analysis/mom/patchy_NA/Analysis/levitus_ewg_temp_salt_cm2m.nc','y_T');

# observation from AVISO product from October 1992 to December 2010
zos = nc_read('/export/grunchfs/unibjerknes/milicak/bckup/obs/zos_AVISO_L4_199210-201012.nc','zos')
lonclim = nc_read('/export/grunchfs/unibjerknes/milicak/bckup/obs/zos_AVISO_L4_199210-201012.nc','lon')
latclim = nc_read('/export/grunchfs/unibjerknes/milicak/bckup/obs/zos_AVISO_L4_199210-201012.nc','lat')
sshclim = np.copy(zos)
sshclim = sshclim[3:-12]
sshclim[sshclim>10] = np.nan
sshclim = np.nanmean(sshclim, axis=0)
map_file = '/bcmhsm/milicak/RUNS/noresm/CORE2/Arctic/DATA/maps/map_woa09_1deg_to_mom_patch_.nc';
n_a = ncgetdim(map_file,'n_a');
n_b = ncgetdim(map_file,'n_b');
# mom 1 degree grid points
nx_b = 360
ny_b = 200
sshclim_src = sshclim.flatten()
mapping_data=io.loadmat('/fimm/home/bjerknes/milicak/Analysis/mom/patchy_NA/Analysis/map_woa09_1deg_to_mom_1deg_patch.mat')
Sspr = mapping_data['Sspr'][:]
sshclim_mom = ma.reshape(Sspr*sshclim_src,(int(ny_b), int(nx_b)))


fig = plt.figure()
#m=Basemap(llcrnrlon=0,llcrnrlat=-88,urcrnrlon=360,urcrnrlat=88,projection='cyl')
m=Basemap(llcrnrlon=-280,llcrnrlat=-88,urcrnrlon=80,urcrnrlat=88,projection='cyl')
m.drawcoastlines()
m.fillcontinents()
m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
#im1 = m.pcolormesh(lonclim,latclim,np.ma.masked_invalid(sshclim),shading='flat',cmap=cmap_needjet2
#                  ,vmin=-1.5,vmax=1.5,latlon=True)
im1 = m.pcolormesh(lon,lat,np.ma.masked_invalid(sshclim_mom),shading='flat',cmap=cmap_needjet2
                  ,vmin=-1.5,vmax=1.5,latlon=True)
cb = m.colorbar(im1,"right", size="5%", pad="10%") # pad is the distance between colorbar and figure
cb.set_label('[m]')
plt.savefig('paperfigs/ssh_clim.eps', bbox_inches='tight',format='eps', dpi=300)
plt.clf()
plt.close(fig)
sys.exit()

for i, project in enumerate(projects):
    filename=root_folder+project+'/'+hist_folder[0]+'/'+'00630101.ocean_month.nc'
    print filename
    #ssh = ncread_time_surface(filename,'eta_t',700,701,nx,ny)
    #ssh = ncread_time_surface(filename,'eta_t',348,744,nx,ny)
    ssh = ncread_time_surface(filename,'eta_t',516,744,nx,ny)
    if i==0:
        ssh_ctrl=np.copy(ssh)
    if i==2:
        ssh2=np.copy(ssh)


    fig = plt.figure()
    m=Basemap(llcrnrlon=-280,llcrnrlat=-88,urcrnrlon=80,urcrnrlat=88,projection='cyl')
    #m = Basemap(projection='merc',llcrnrlat=-80,urcrnrlat=80,\
    #        llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c')
    m.drawcoastlines()
    m.fillcontinents()
    m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
    m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
    im1 = m.pcolormesh(lon,lat,np.ma.masked_invalid(ssh+0.5),shading='flat',cmap=cmap_needjet2
                      ,vmin=-1.5,vmax=1.5,latlon=True)
    cb = m.colorbar(im1,"right", size="5%", pad="10%") # pad is the distance between colorbar and figure
    cb.set_label('[m]')
    plt.show()
    plt.savefig('paperfigs/'+project+'_ssh.eps', bbox_inches='tight',format='eps', dpi=300)
    plt.clf()
    plt.close(fig)
    fig = plt.figure()
    m=Basemap(llcrnrlon=-280,llcrnrlat=-88,urcrnrlon=80,urcrnrlat=88,projection='cyl')
    #m = Basemap(projection='merc',llcrnrlat=-80,urcrnrlat=80,\
    #        llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c')
    m.drawcoastlines()
    m.fillcontinents()
    m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
    m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
    im1 = m.pcolormesh(lon,lat,np.ma.masked_invalid(ssh+0.5-sshclim_mom),
                       shading='flat',cmap=discrete_cmap(32, 'RdBu_r')
                      ,vmin=-.5,vmax=.5,latlon=True)
    cb = m.colorbar(im1,"right", size="5%", pad="10%") # pad is the distance between colorbar and figure
    cb.set_label('[m]')
    plt.show()
    plt.savefig('paperfigs/'+project+'_sshbias.eps', bbox_inches='tight',format='eps', dpi=300)
    plt.clf()
    plt.close(fig)



# difference from sst_ctrl
fig = plt.figure()
m=Basemap(llcrnrlon=-280,llcrnrlat=-88,urcrnrlon=80,urcrnrlat=88,projection='cyl')
# llcrnrlat,llcrnrlon,urcrnrlat,urcrnrlon
# are the lat/lon values of the lower left and upper right corners
# of the map.
# lat_ts is the latitude of true scale.
# resolution = 'c' means use crude resolution coastlines.
#m = Basemap(projection='merc',llcrnrlat=-80,urcrnrlat=80,\
#            llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c')
m.drawcoastlines()
m.fillcontinents()
m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
im1 = m.pcolormesh(lon,lat,np.ma.masked_invalid(ssh2-ssh_ctrl),shading='flat',cmap=discrete_cmap(32, 'RdBu_r')
                    ,vmin=-.25,vmax=.25,latlon=True)
cb = m.colorbar(im1,"right", size="5%", pad="10%") # pad is the distance between colorbar and figure
cb.set_label('[m]')
plt.show()
plt.savefig('paperfigs/'+project+'_ssh_diff_patchy2_ctrl.eps', bbox_inches='tight',format='eps', dpi=300)
plt.clf()
plt.close(fig)
