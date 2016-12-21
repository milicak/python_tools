import numpy as np
#%matplotlib inline
#np.shape !!!!!
from mpl_toolkits.basemap import Basemap
# another way to plot
#import cartopy.crs as ccrs
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
from eofs.standard import Eof

grid_file='/mnt/fimm/Analysis/pop/Arctic/Analysis/areacello_fx_CCSM4_piControl_r0i0p0.nc';


root_folder='/mnt/fimm/Analysis/pop/Arctic/Analysis/matfiles/'

projects=['B1850CN_f19_tn11_kdsens','B1850CN_f19_tn11_kdsens01','B1850CN_f19_tn11_kdsens02','B1850CN_f19_tn11_kdsens03',
          'B1850CN_f19_tn11_kdsens05','B1850CN_f19_tn11_kdsens04','B1850CN_f19_tn11_kdsens06']

#05 ctrl for warm; 04 Barents,Kara + Eurasia; 06 Barents,Kara for warm exps
nx = 144;
ny = 96;

av1=75; #240; %75;
av2=85; #250; %85;

fyear = '1'; # first year
lyear = '250'; # last year
south1 = 20; #20 North
north1 = 80; # 80 North
west1 = -80+360; #90 West
east1 = 40; #40 East

for i in xrange(0,1):                 
    print "mehmet1",i,projects[i]
    filename=root_folder+projects[i]+'_airslp_mean_'+fyear+'_'+lyear+'.mat'
    mat = scipy.io.loadmat(filename)
    lon=np.array(mat['lon'])
    lat=np.array(mat['lat']) 
    [lons,lats]=np.meshgrid(lon,lat)
    coslat = np.cos(np.deg2rad(lat)).clip(0., 1.)
    wgts = np.sqrt(coslat)[..., np.newaxis]
    jind1=np.max(np.where(lat<=south1))
    jind2=np.max(np.where(lat<=north1))
    iind1=np.max(np.where(lon<=west1))
    iind2=np.max(np.where(lon<=east1))
    mask=np.zeros([np.int(lyear),nx,ny])*np.nan
    mask[:,iind1:,jind1:jind2]=1
    mask[:,0:iind2-1,jind1:jind2]=1
    slp_winter=np.array(mat['AIRslpmeanwinter'])
    data=slp_winter*mask
    # remove mean
    data_mean=np.nanmean(data,0)
    data=data-data_mean;
    # Compute anomalies by removing the time-mean.
    #data_mean = data.mean(axis=0)
    #data=data-data_mean;
    solver = Eof(data, weights=wgts)

    # Retrieve the leading EOF, expressed as the covariance between the leading PC
    # time series and the input SLP anomalies at each grid point.
    eofS = solver.eofsAsCovariance(neofs=3)
    pcS = solver.pcs(npcs=3, pcscaling=1)
    variance_eofs=solver.varianceFraction(neigs=3)
    #m = Basemap(width=8000000,height=8000000,
    #        resolution='l',projection='stere',\
    #        lat_ts=40,lat_0=90,lon_0=0.)
    m = Basemap(projection='cyl',llcrnrlat=0,urcrnrlat=90,\
            llcrnrlon=-100,urcrnrlon=10,resolution='c')
                        
    m.drawcoastlines() #m.drawcoastlines(linewidth=1.5)
    m.fillcontinents()
    #m.drawparallels(parallels)
    #m.drawmeridians(meridians)
    # Plot the leading EOF expressed as covariance in the European/Atlantic domain.
    im1 = m.pcolormesh(lons,lats,np.transpose(np.ma.masked_invalid(np.squeeze(eofS[2,:,:]))),shading='flat',cmap='jet',latlon=True)
    #plt.title('EOF1 expressed as covariance', fontsize=16)
    plt.show()