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

# sdate = "%c%4.4d%c" % ('*',fyear,'*')
fname = root_folder+project_name+'/'+expid+'/OUT/'+'uTSS_lobc_chunk_'+'*'+'.ous.nc'
grd = xr.open_dataset(root_folder+project_name+'/'+expid+'/OUT/'+'uTSS_lobc_chunk_'+'0001'+'.nos.nc')
grd['element_index'] -= 1

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


# shapef,area = compute_shapefnc(grd)

list = sorted(glob.glob(fname))
# dfu = xr.open_mfdataset(list)['u_velocity'][:,:,0] # SST
# dfv = xr.open_mfdataset(list)['v_velocity'][:,:,0] # SST
ind1 = 59 # 100 m to 500 m
ind2 = 71 # 100 m to 500 m

ind1 = 40 # 40 m to 100 m
ind2 = 59 # 40 m to 100 m
dfu = xr.open_mfdataset(list)['u_velocity'][:,:,ind1:ind2] 
dfv = xr.open_mfdataset(list)['v_velocity'][:,:,ind1:ind2] # 100 m to 500 m
zw = 0.5*(np.copy(grd.level[ind1:ind2+1])+np.copy(grd.level[ind1-1:ind2]))
dz = zw[1:]-zw[:-1]
# dfu = xr.open_mfdataset(list)['u_velocity'][:,:,40:59] # SST
# dfv = xr.open_mfdataset(list)['v_velocity'][:,:,40:59] # SST
## # dfu = dfu.mean('time')
## # dfv = dfv.mean('time')
sp = np.sqrt(dfu**2+dfv**2)
sp = sp.mean('time')
mask1 = np.ones([sp.shape[0],sp.shape[1]])
mask1[sp==0] = 0
dnm = sp*mask1*dz
dnm2 = mask1*dz
spmean = np.nansum(dnm,1)/np.nansum(dnm2,1)

# dfun=(1.0/3.0)*(dfu[grd.element_index[:,0]] + dfu[grd.element_index[:,1]] + dfu[grd.element_index[:,2]])
# dfvn=(1.0/3.0)*(dfv[grd.element_index[:,0]] + dfv[grd.element_index[:,1]] + dfv[grd.element_index[:,2]])


# df['time'] = pd.date_range('2016-01-01', periods=318)
# df_season = df.where(df['time.season'] == 'DJF').groupby('time.year').mean('time')


plt.figure()
# m = Basemap(llcrnrlon=22.5,llcrnrlat=38.5,urcrnrlon=32.,urcrnrlat=43.5,\
m = Basemap(llcrnrlon=26.0,llcrnrlat=40.0,urcrnrlon=30.,urcrnrlat=41.4,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='h',projection='merc',\
            lat_0=40.,lon_0=20.,lat_ts=20.)

m.drawcoastlines(linewidth=0.2)
m.fillcontinents(color='grey')
m.drawparallels(np.arange(38,44,1),labels=[1,1,0,0])
m.drawmeridians(np.arange(22,33,2),labels=[0,0,0,1])
longitude,latitude = m(np.copy(grd.longitude),np.copy(grd.latitude))

point_mask=~np.isfinite(spmean)
# triang = mtri.Triangulation(longitude,latitude)
# tri_mask = np.any(point_mask.data[triang.triangles],axis=1)
tri_mask = np.any(point_mask[grd.element_index],axis=1)

# im1=plt.tripcolor(longitude,latitude,grd.element_index,
                # sp,cmap='needJet2',vmin=0,vmax=0.6,shading='gouraud')
im1=plt.tripcolor(longitude,latitude,grd.element_index,
                spmean,cmap='needJet2',mask=tri_mask,vmin=0,vmax=0.05,shading='gouraud')
cb = m.colorbar(im1,"right", size="5%", pad="10%", extend='max')



