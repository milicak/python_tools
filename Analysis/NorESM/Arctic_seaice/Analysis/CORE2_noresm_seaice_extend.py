#!/usr/bin/env python
import numpy as np
import numpy.ma as ma
import pdb
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import sys
from matplotlib.path import Path as mpPath
from netcdf_functions import nc_read
from netcdf_functions import ncgetdim
import NorESM_utils as noresmutils
import scipy.io

plt.ion()

datesep = '-'


# global constants
aradius = 6.373E6      # Radius of Earth (m)
parallels = np.arange(-80.,90,20.)
meridians = np.arange(0.,360.,20.)

def enable_global(tlon,tlat):
  """Fix the data in such a way that it can to be plotted on a global projection on its native grid"""
  tlon = np.where(np.greater_equal(tlon,min(tlon[:,0])),tlon-360,tlon)
  tlon=tlon+abs(ma.max(tlon)); tlon=tlon+360
  # stack grids side-by-side (in longitiudinal direction), so
  # any range of longitudes may be plotted on a world map.
  tlon = np.concatenate((tlon,tlon+360),1)
  tlat = np.concatenate((tlat,tlat),1)
  tlon = tlon-360.
  return tlon, tlat

def ncread_time_surface(fname, variable, timestr, timeend, x, y):
    # how to use this subroutine is from netcdf_functions import nc_read
    ncfile = Dataset(fname, 'r', format='NETCDF4')
    tmp = np.zeros([y, x])
    for i in range(timestr, timeend):
        # print i
        tmp = tmp+ncfile.variables[variable][i, :, :].copy()

    tmp = tmp/(timeend-timestr)
    return tmp


def ncread_lev(fname, variable, zlev):
    ncfile = Dataset(fname, 'r', format='NETCDF4')
    tmp = np.squeeze(ncfile.variables[variable][:,zlev,:,:].copy())

    return tmp


def get_sdate_ini(root_folder,cmpnt,mdl,ext, **kwargs):
    expid = kwargs.get('expid', None)
    m2y = kwargs.get('m2y', None)
    # Get dimensions and time attributes for ocn or atm or others ...
    prefix=root_folder+expid+'/'+cmpnt+'/hist/'+expid+ '.'+ mdl + \
            '.'+ext+'.'
    if m2y==1:
        sdate="%4.4d%c%2.2d" % (1,'-',1)    # assuming fyear 1 exists
    else:
        sdate="%4.4d" % (1)


    return prefix, sdate


def timemean(prefix, var, varname, fyear, lyear, **kwargs):
    zlev = kwargs.get('zlev', None)
    expid = kwargs.get('expid', None)
    m2y = kwargs.get('m2y', None)
    ''' timemean averages '''
    months2days=[31,  28,  31,  30,  31,   30,   31,  31,   30, 31,   30, 31];
    yeardays=sum(months2days);
    mw = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31], dtype=np.float)
    mw = mw/sum(mw)
    n=0.0;
    fyear
    lyear
    for year in xrange(fyear,lyear+1):
        if m2y==1:
            for month in xrange(0,12):
                n=n+mw[month]
                sdate="%4.4d%c%2.2d" % (year,'-',month+1)
                if not zlev == None:
                    dnm=ncread_lev(prefix+sdate+'.nc', varname, zlev);
                else:
                    dnm=nc_read(prefix+sdate+'.nc', varname);


                var=var+np.squeeze(dnm.data)*mw[month]
            print sdate, n


        else:
            sdate = "%4.4d" % (year)
            n += 1.0
            if not zlev == None:
                dnm=ncread_lev(prefix+sdate+'.nc', varname, zlev);
            else:
                dnm=nc_read(prefix+sdate+'.nc', varname);


            var=var+np.squeeze(dnm.data)
            print sdate, n


    var = var/n
    return var


def main():
    # tripolar 1degree grid
    filename = "/fimm/home/bjerknes/milicak/Analysis/NorESM/Arctic_seaice/Analysis/region_masks.mat"
    mat = scipy.io.loadmat(filename)
    # lon1,lat1 is for Kara and Barents Sea
    # lon2,lat2 is for Greenland Sea
    # lon3,lat3 is for Hudson Bay
    # lon4,lat4 is for CAA
    # lon5,lat5 is for Arctic Ocean Canadian side
    # lon6,lat6 is for Labrador Sea/ Baffin Bay
    # lon7,lat7 is for Arctic Ocean Eurasian side
    # lon8,lat8 is for Bering Sea
    lonv = np.array(mat['lon8'])
    latv = np.array(mat['lat8'])

    grid_file='/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid.nc';
    # tripolar 0.25degree grid
    #grid_file = '/bcmhsm/milicak/RUNS/noresm/CORE2/Arctic/maps/grid_0_25degree.nc';
    # bi-polar grid
    #grid_file='/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid_bipolar.nc';
    plt.figure()
    m = Basemap(width=8000000,height=8000000,
            resolution='l',projection='stere',\
            lat_ts=40,lat_0=90,lon_0=0.)
    lon = nc_read(grid_file,'plon')
    lat = nc_read(grid_file,'plat')
    # doesn't work for lon2
    # lon2 no change!
    lon[lon<-10] = lon[lon<-10]+360
    lonv[lonv<-10] = lonv[lonv<-10]+360
    # for lon 7
    #lon[lon<-110] = lon[lon<-110]+360
    #lonv[lonv<-110] = lonv[lonv<-110]+360
    #[lon,lat] = enable_global(lon,lat)
    nx = lon.shape[1]
    ny = lon.shape[0]
    xpt,ypt = m(lonv,latv)
    vertices = np.transpose(np.array([lonv.flatten(),latv.flatten()]))
    # vertices = [(146.0,-42.0),(167.0,-42.0),(167.0,-53.0), (146.0,-53.0)]
    pnts = np.transpose(np.array([lon.flatten(),lat.flatten()]))
    path = mpPath(vertices)
    inside = path.contains_points(pnts)
    inside = np.reshape(inside,((ny,nx)))
    m.drawcoastlines() #m.drawcoastlines(linewidth=1.5)
    #m.fillcontinents()
    m.drawparallels(parallels)
    m.drawmeridians(meridians)
    im1 = m.pcolormesh(lon,lat,np.ma.masked_invalid(inside),
                       shading='flat',latlon=True)
    m.plot(xpt,ypt,color='r')
    m.plot(xpt[0],ypt[0],'ko')
    m.plot(xpt[1],ypt[1],'go')
    m.plot(xpt[2],ypt[2],'mo')
    m.plot(xpt[3],ypt[3],'yo')
    m.plot(xpt[4],ypt[4],'ko')
    m.plot(xpt[5],ypt[5],'co')
    print 'mehmet',lonv
    print 'mehmet2',latv





if __name__ == "__main__":
    # this won't be run when imported
    main()

