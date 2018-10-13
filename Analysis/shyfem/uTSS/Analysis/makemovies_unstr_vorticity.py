import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib import rc
# import ESMF
from mpl_toolkits.basemap import Basemap                                            
#import cartopy.crs as ccrs                                                          
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
 
# plt.ion()
plt.ioff()
#rc('text', usetex=True)

root_folder  = '/work/mi19918/Projects/'
project_name = 'uTSS'

expid = 'Exp_20160101'
#expid = 'Exp01.2'

Nnode = 240993
# coriolis(41N)
fzero = 9.5681e-5

fyear = 0;
lyear = 30;

sdate = "%c%4.4d%c" % ('*',fyear,'*')
fname = root_folder+project_name+'/'+expid+'/OUT/'+'uTSS_lobc_chunk_'+sdate+'.ous.nc'
list = sorted(glob.glob(fname))
str1 = ''.join(list)
grd = xr.open_dataset(str1)
for year in xrange(fyear+1,lyear+1):
    sdate = "%c%4.4d%c" % ('*',year,'*')
    list.extend(sorted(glob.glob(root_folder+project_name+'/'+expid+'/OUT/'+'uTSS_lobc_chunk_'+sdate+'.ous.nc')))



data = xr.open_mfdataset(list, chunks={'time':5, 'node':2000})
data['element_index'] -= 1
grd['element_index'] -= 1
#data = xr.open_dataset(fname)['salinity']


# triang = mtri.Triangulation(data.longitude,data.latitude,data.element_index)


def coriolis(lat):  
    """Compute the Coriolis parameter for the given latitude:
    ``f = 2*omega*sin(lat)``, where omega is the angular velocity 
    of the Earth.
    
    Parameters
    ----------
    lat : array
      Latitude [degrees].
    """
    omega   = 7.2921159e-05  # angular velocity of the Earth [rad/s]
    return 2*omega*np.sin(lat/360.*2*np.pi)


def ev_g2c(lambdar, phi, lambda0, phi0):
    ''' dadsasd   '''
    rad = np.pi/180.0
    rEarth = 6378206.4E0
    dlambda = rad * (lambdar - lambda0)                               
    dphi    = rad * (phi - phi0)                                     
    xa = rEarth*dlambda*np.cos(rad*phi0)                                     
    ya = rEarth*dphi      

    return xa,ya



shapef = {}
def compute_shapefnc(grd):
    dlon0 = 26.68237 
    dlat0 = 40.54126
    x1tmp = grd.longitude[grd.element_index[:,0]]
    y1tmp = grd.latitude[grd.element_index[:,0]]
    x1,y1 = ev_g2c(x1tmp, y1tmp, dlon0, dlat0)
    x2tmp = grd.longitude[grd.element_index[:,1]]
    y2tmp = grd.latitude[grd.element_index[:,1]]
    x2,y2 = ev_g2c(x2tmp, y2tmp, dlon0, dlat0)
    x3tmp = grd.longitude[grd.element_index[:,2]]
    y3tmp = grd.latitude[grd.element_index[:,2]]
    x3,y3 = ev_g2c(x3tmp, y3tmp, dlon0, dlat0)
    area = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
    aji = 1.0/area
    b1 = (y2-y3)*aji
    c1 = (x3-x2)*aji
    b2 = (y3-y1)*aji
    c2 = (x1-x3)*aji
    b3 = (y1-y2)*aji
    c3 = (x2-x1)*aji
    dnm = [b1,b2,b3]
    shapef.setdefault('xderiv',[]).append(dnm)
    dnm = [c1,c2,c3]
    shapef.setdefault('yderiv',[]).append(dnm)
    return shapef,area


shapef,area = compute_shapefnc(grd)

# velocity time,node,level indexing
varu = data.u_velocity[:,:,0]
varv = data.v_velocity[:,:,0]

u1 = varu[:,grd.element_index[:,0]]
u2 = varu[:,grd.element_index[:,1]]
u3 = varu[:,grd.element_index[:,2]]
v1 = varv[:,grd.element_index[:,0]]
v2 = varv[:,grd.element_index[:,1]]
v3 = varv[:,grd.element_index[:,2]]

# relative vorticity ksi = vx-uy
vx = v1*shapef['xderiv'][0][0]
vx = vx+v2*shapef['xderiv'][0][1]
vx = vx+v3*shapef['xderiv'][0][2]

uy = u1*shapef['yderiv'][0][0]
uy = uy+u2*shapef['yderiv'][0][1]
uy = uy+u3*shapef['yderiv'][0][2]

ksi = vx-uy

# compute ksi at node points
ksinode=np.zeros((ksi.shape[0],Nnode))
ksinode[:,grd.element_index[grd.element,0]] = ksinode[:,grd.element_index[grd.element,0]]+ksi[:,grd.element]*area
ksinode[:,grd.element_index[grd.element,1]] = ksinode[:,grd.element_index[grd.element,1]]+ksi[:,grd.element]*area
ksinode[:,grd.element_index[grd.element,2]] = ksinode[:,grd.element_index[grd.element,2]]+ksi[:,grd.element]*area

weight = np.zeros(Nnode)
weight[grd.element_index[grd.element,0]] = weight[grd.element_index[grd.element,0]]+area
weight[grd.element_index[grd.element,1]] = weight[grd.element_index[grd.element,1]]+area
weight[grd.element_index[grd.element,2]] = weight[grd.element_index[grd.element,2]]+area

ksinode = ksinode/weight

fcor = coriolis(grd.latitude)

for ind in xrange(fyear,lyear+1):
    print ind
    sdate = "%4.4d" % (ind)
    plt.figure(figsize=(12,8))
    m = Basemap(llcrnrlon=22.5,llcrnrlat=38.5,urcrnrlon=32.,urcrnrlat=43.5,\
                rsphere=(6378137.00,6356752.3142),\
                resolution='h',projection='merc',\
                lat_0=40.,lon_0=20.,lat_ts=20.)
    
    m.drawcoastlines(linewidth=0.2)
    m.fillcontinents(color='grey')
    m.drawparallels(np.arange(38,44,1),labels=[1,1,0,0])
    m.drawmeridians(np.arange(22,33,2),labels=[0,0,0,1])
    
    longitude,latitude = m(np.copy(grd.longitude),np.copy(grd.latitude))
    im1=plt.tripcolor(longitude,latitude,grd.element_index,
                      ksinode[ind,:]/fcor,cmap='seismic',vmin=-1,vmax=1,shading='gouraud')
    cb = m.colorbar(im1,"right", size="5%", pad="10%",ticks=[-1,-0.75,-0.5, -0.25,
                                                             0, 0.25, 0.5, 0.75, 1]) # pad is the distance between colorbar and figure
    cb.set_label('$\zeta/f$',rotation=0,y=1.07,labelpad=-45)
    printname = 'gifs/vorticity_'+sdate+'.png'
    plt.savefig(printname, bbox_inches='tight',format='png',dpi=300)
    plt.close()




