import my_nanfilter
import numpy as np
import numpy.ma as ma
import scipy.io
import sys
import os
#%matplotlib inline
#np.shape !!!!!
from needJet2 import shfn
from mpl_toolkits.basemap import Basemap
from cpttoseg import cpt2seg
from my_nanfilter import my_nanfilterbox
from netCDF4 import Dataset
from netcdf_functions import nc_read
from netcdf_functions import ncgetdim
import matplotlib.pyplot as plt
import matplotlib as mpllib
reload(my_nanfilter)

from scipy.sparse import csr_matrix

# IMPORTANT 
plt.ion()

lon=nc_read('/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid.nc','plon')
lat=nc_read('/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid.nc','plat')
lon=np.copy(lon[:-1,:])
lat=np.copy(lat[:-1,:])

cmap_needjet2=shfn()

_palette_data = cpt2seg('/fimm/home/bjerknes/milicak/python_tools/Analysis/cpt_files/sst.cpt')
palette = mpllib.colors.LinearSegmentedColormap('palette', _palette_data, 256)

def ncread_time(fname,variable,timestr,timeend,x,y):
    # how to use this subroutine is from netcdf_functions import nc_read
    ncfile = Dataset(fname, 'r', format='NETCDF4')
    tmp=np.zeros([y,x])
    for i in range(timestr,timeend):
        #print i
        tmp=tmp+ncfile.variables[variable][i,:,:].copy();

    tmp=tmp/(timeend-timestr)
    return tmp


def enable_global(tlon,tlat,data):
  """Fix the data in such a way that it can to be plotted on a global projection on its native grid"""
  tlon = np.where(np.greater_equal(tlon,min(tlon[:,0])),tlon-360,tlon)
  tlon=tlon+abs(ma.max(tlon)); tlon=tlon+360
  # stack grids side-by-side (in longitiudinal direction), so
  # any range of longitudes may be plotted on a world map.
  tlon = np.concatenate((tlon,tlon+360),1)
  tlat = np.concatenate((tlat,tlat),1)
  data = ma.concatenate((data,data),1)
  tlon = tlon-360.
  return tlon, tlat, data


fyear = '650'; # first year
lyear = '750'; # last year
nx = 360;
ny = 385;
nz = 53;

i0=25;
j0=290;

root_folder='/work/milicak/mnt/norstore/NS2345K/noresm/cases/'

projects=['N1850_f19_tn11_reverseCoriolis']
#projects=['N1850_f19_tn11_01_default']

uflx=np.zeros([ny,nx])
vflx=np.zeros([ny,nx])
    
for year in xrange(np.int(fyear),np.int(lyear)+1):
    #for month in xrange(1,13):
    #filename = root_folder+projects[0]+'atm/hist/'+projects[0]+'cam.h0.'+str(year).zfill(4)+'-'+str(month).zfill(2)+'.nc'
    filename = root_folder+projects[0]+'/ocn/hist/'+projects[0]+'.micom.hy.'+str(year).zfill(4)+'.nc'
    dnmu = np.squeeze(nc_read(filename,'uflx')*1e-9)
    dnmv = np.squeeze(nc_read(filename,'vflx')*1e-9)
    dnmu = np.squeeze(np.sum(dnmu,0))
    dnmv = np.squeeze(np.sum(dnmv,0))
    uflx = uflx + dnmu 
    vflx = vflx + dnmv
    print year

uflx = uflx/(np.float(lyear)-np.float(fyear)+1)
vflx = vflx/(np.float(lyear)-np.float(fyear)+1)


# compute stream function following Mats Bentsen ;

qi=np.linspace(1,nx*ny,nx*ny);
qis=np.concatenate((np.arange((nx+1),2*nx+1,1),np.arange(1,nx*(ny-1)+1,1)));
qin=np.concatenate((np.arange((nx+1),(nx*(ny-1)+nx/2+1)+1,1)
                   ,np.arange((nx*(ny-1)+nx/2),(nx*(ny-1)+2)-1,-1)));
dnm=np.concatenate(((nx*(ny-1)+1,),np.arange(nx*(ny-1),(nx*(ny-2)+nx/2+2)-1,-1)))                 
qin=np.concatenate((qin,dnm));
dnm=np.concatenate((((nx*(ny-1)+nx/2+1),),np.arange((nx*(ny-2)+nx/2),(nx*(ny-2)+2)-1,-1)))
qin=np.concatenate((qin,dnm));

x=np.reshape(qi,(ny,nx))
xx=np.roll(x,1,axis=1);
qie=np.reshape(xx,(nx*ny,))

x=np.reshape(qi,(ny,nx))
xx=np.roll(x,-1,axis=1);
qiw=np.reshape(xx,(nx*ny,))

wi=np.arange((nx+1),nx*ny+1,1);
uis=np.arange(1,nx*(ny-1)+1,1);
x=np.reshape(wi,(ny-1,nx))
xx=np.roll(x,1,axis=1);
vie=np.reshape(xx,(nx*(ny-1),))

row=np.concatenate((qi,qi,qi,qi,qi));
col=np.concatenate((qi,qis,qie,qiw,qin));
x4=4*np.ones(nx*ny,)
x1=-np.ones(4*nx*ny,)
s=np.concatenate((x4,x1))

omega=np.zeros(nx*ny,);
omega[wi-1]=-np.reshape(uflx,(nx*ny)).data[uis-1]+np.reshape(uflx,(nx*ny)).data[wi-1]+np.reshape(vflx,(nx*ny)).data[vie-1]-np.reshape(vflx,(nx*ny)).data[wi-1]
#sys.exit()

# creating A matrix is not necessary because np.linalg.solve crashes 
# So we will save row, col, s and omega to txt and call a matlab function to
# solve streamfunction equation.

np.savetxt('strm_s.txt',s)
np.savetxt('strm_row.txt',row)
np.savetxt('strm_col.txt',col)
np.savetxt('strm_omega.txt',omega)

#sys.exit()
print 'calling matlab now'
os.system('/fimm/apps/matlab/2013a/bin/matlab -nodesktop -nojvm < compute_strm_micom.m')


mat = scipy.io.loadmat('strm_strmf.mat')
strmf=np.array(mat['strmf'])
strmf=np.copy(strmf[:,:-1])
[lon,lat,strmf]=enable_global(lon,lat,np.transpose(strmf))

#A=csr_matrix((s,(row-1, col-1)), shape=(nx*ny, nx*ny)).toarray();
#A=scipy.sparse.csr_matrix((s,(row-1, col-1)), shape=(nx*ny, nx*ny),dtype=np.int8);
#A=scipy.sparse.csr_matrix((s,(row-1, col-1)), shape=(nx*ny,
#coeff_mat = scipy.sparse.coo_matrix((s, (row-1, col-1))).tocsc()
#dnm=np.linalg.solve(A,omega)

fig = plt.figure()
#ax = plt.gca()
#ax.set_axis_bgcolor('grey')
m=Basemap(llcrnrlon=-180,llcrnrlat=-88,urcrnrlon=180,urcrnrlat=88,projection='cyl')
m.drawcoastlines()
m.fillcontinents()
m.drawparallels(np.arange(-80,81,20),labels=[1,1,0,0])
m.drawmeridians(np.arange(0,360,60),labels=[0,0,0,1])
im1 =m.pcolormesh(np.transpose(lon-110),np.transpose(lat),
                  np.transpose(np.ma.masked_invalid(strmf))
                  ,shading='flat',cmap=cmap_needjet2,vmin=-200,vmax=70,latlon=True)
                  #,shading='flat',cmap=cmap_needjet2,vmin=-150,vmax=150,latlon=True)
cb = m.colorbar(im1,"right", size="5%", pad="10%") #,ticks=[-4, -3, -2, -1, 0, 1, 2, 3, 4]) #dnm pad is the distance between colorbar and figure
cb.set_label('[Sv]')
#plt.ylabel('Depth [m]')
#plt.xlabel('Lat')
#plt.ylim(-7000,0)
#cb.set_label('[' r'$^\circ$' 'C]')
#sys.exit()
plt.savefig('paperfigs/'+projects[0]+'_streamfunction_mean.eps', 
            bbox_inches='tight',format='eps', dpi=200)
#plt.clf()
#plt.close(fig)


#plt.show()



