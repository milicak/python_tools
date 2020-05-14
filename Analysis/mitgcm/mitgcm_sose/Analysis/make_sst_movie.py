import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import xmitgcm
import os
from datetime import date
import matplotlib.pyplot as plt
# import ESMF
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LinearSegmentedColormap
import math
plt.ioff()

root_folder = '/cluster/work/users/milicak/RUNS/mitgcm/mitgcm_sose/'

expid = 'Exp01_0';

sshname = 'Southern_Ocean_ctrl_daily_ocn_ETAN_2009_03_01cyc.nc'
sstname = 'Southern_Ocean_ctrl_daily_ocn_SST_2009_03_01cyc.nc'
vortname = 'Southern_Ocean_ctrl_daily_ocn_vort2dsurf_2009_03_01cyc.nc'

fname = root_folder+expid+'/'+sstname
dt = xr.open_dataset(fname)

fname = root_folder+expid+'/'+vortname
dv = xr.open_dataset(fname)

fname = root_folder+expid+'/'+sshname
ds = xr.open_dataset(fname)

lonc, latc = np.meshgrid(np.copy(dt.XC),np.copy(dt.YC))

def polar_stere(lon_w, lon_e, lat_s, lat_n, **kwargs):
    '''Returns a Basemap object (NPS/SPS) focused in a region.

    lon_w, lon_e, lat_s, lat_n -- Graphic limits in geographical coordinates.
                                  W and S directions are negative.
    **kwargs -- Aditional arguments for Basemap object.

    '''
    lon_0 = lon_w + (lon_e - lon_w) / 2.
    ref = lat_s if abs(lat_s) > abs(lat_n) else lat_n
    lat_0 = math.copysign(90., ref)
    proj = 'npstere' if lat_0 > 0 else 'spstere'
    prj = Basemap(projection=proj, lon_0=lon_0, lat_0=lat_0,
                          boundinglat=0, resolution='c')
    #prj = pyproj.Proj(proj='stere', lon_0=lon_0, lat_0=lat_0)
    lons = [lon_w, lon_e, lon_w, lon_e, lon_0, lon_0]
    lats = [lat_s, lat_s, lat_n, lat_n, lat_s, lat_n]
    x, y = prj(lons, lats)
    ll_lon, ll_lat = prj(min(x), min(y), inverse=True)
    ur_lon, ur_lat = prj(max(x), max(y), inverse=True)
    return Basemap(projection='stere', lat_0=lat_0, lon_0=lon_0,
                           llcrnrlon=ll_lon, llcrnrlat=ll_lat,
                           urcrnrlon=ur_lon, urcrnrlat=ur_lat, **kwargs)



for ind in range(0,30):
    sdate="%2.2d" % (ind+1)
    print(sdate)
    plt.figure(figsize=(10,10))
    m = Basemap(projection='spstere',boundinglat=-25,lon_0=0,resolution='l')
    # lonstr = -40+ind
    # m= polar_stere(lonstr, lonstr+50.0, -88, -25, resolution='l')
    # m= polar_stere(-5.0, 45.0, -88, -25, resolution='l')
    lon, lat = m(lonc,latc);
    m.bluemarble();
    # m.pcolormesh(lon,lat,ma.masked_equal(dt.THETA[0,:,:],0),cmap='nice_gfdl',
                 # vmin=-2,vmax=28);
    #
    m.pcolormesh(lon[:,:2160],lat[:,:2160],
                 ma.masked_equal(dt.THETA[ind,:,:2160],0),
                 cmap='nice_gfdl',vmin=-2,vmax=28);
    m.pcolormesh(lon[:,2160:],lat[:,2160:],
                 ma.masked_equal(dv.momVort3[ind,:,2160:],0),
                 cmap='RdBu_r',vmin=-1e-4,vmax=1e-4);
    savename = 'gifs/sst_vort_'+sdate+'.png'
    plt.savefig(savename,bbox_inches='tight',format='png',dpi=150)
    plt.clf()
    plt.close()

    # m.drawcoastlines()
    # m.fillcontinents(color='coral')
    # m.pcolormesh(lon,lat,ma.masked_equal(ds.SIarea,0),cmap=cmocean.cm.ice);plt.colorbar()
