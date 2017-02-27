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

plt.ion()

datesep = '-'


# global constants
aradius = 6.373E6      # Radius of Earth (m)

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

    grid_file='/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid.nc';
    # tripolar 0.25degree grid
    #grid_file = '/bcmhsm/milicak/RUNS/noresm/CORE2/Arctic/maps/grid_0_25degree.nc';
    # bi-polar grid
    #grid_file='/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid_bipolar.nc';
    lon = nc_read(grid_file,'plon')
    lat = nc_read(grid_file,'plat')
    print 'mehmet',lon.shape
    vertices = [(146.0,-42.0),(167.0,-42.0),(167.0,-53.0), (146.0,-53.0)]
    pnts = np.transpose(np.array([lon.flatten(),lat.flatten()]))
    path = mpPath(vertices)
    inside = path.contains_points(pnts)
    inside = np.reshape(inside,((385,360)))


    #root_folder = '/work/milicak/mnt/norstore/NS2345K/noresm/cases/'
    #root_folder = '/work/milicak/mnt/viljework/archive/'
    root_folder = '/hexagon/work/milicak/archive/'

    #expid = 'NOIIA_T62_tn11_norems2_ctrl_tke'
    #expid = 'NBF1850_f19_tn11_sst_sss_rlx_01'

    gridtype = 'tnx1v1' # tnx1v1 , tnx0.25v1, gx1v6
    woafnamet = '/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/t00an1.nc'
    woafnames = '/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/s00an1.nc'
    #mask_woa09_file='/fimm/home/bjerknes/milicak/Analysis/NorESM/general/Analysis/woa_mask.mat';



if __name__ == "__main__":
    # this won't be run when imported
    main()

