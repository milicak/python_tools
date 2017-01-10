#!/usr/bin/env python 
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import sys
from netcdf_functions import nc_read
from netcdf_functions import ncgetdim

root_folder = '/work/milicak/mnt/norstore/NS2345K/noresm/cases/'

#expid = 'NOIIA_T62_tn11_norems2_ctrl_tke'
expid = 'N1850_f19_tn11_01_default'
fyear = 1; # first year
lyear = 10; # last year

m2y = 0
datesep = '-'

# tripolar 1degree grid
grid_file='/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid.nc';
# tripolar 0.25degree grid
#grid_file = '/bcmhsm/milicak/RUNS/noresm/CORE2/Arctic/maps/grid_0_25degree.nc';
# bi-polar grid
#grid_file='/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid_bipolar.nc';

def get_sdate_ini(args):
    # Get dimensions and time attributes
    if m2y==1:
        prefix=root_folder+expid+'/ocn/hist/'+expid+ '.micom.hm.'
        sdate="%4.4d%c%2.2d" % (fyear,datesep,1)
    else:
        prefix=root_folder+expid+'/ocn/hist/'+expid+ '.micom.hy.'
        sdate="%4.4d" % (fyear)


    return prefix, sdate


def timemean(args):
    ''' timemean averages '''
    mask=nc_read(grid_file,'pmask');

    prefix,sdate = get_sdate_ini(root_folder) 
    months2days=[31,  28,  31,  30,  31,   30,   31,  31,   30, 31,   30, 31];
    yeardays=sum(months2days);

    nx=ncgetdim(prefix+sdate+'.nc','x');
    ny=ncgetdim(prefix+sdate+'.nc','y');
    nz=ncgetdim(prefix+sdate+'.nc','depth');
    depth=nc_read(prefix+sdate+'.nc','depth');

    templvl=np.zeros((nz,ny,nx));
    salnlvl=np.zeros((nz,ny,nx));
    difisolvl=np.zeros((nz,ny,nx));
    difdialvl=np.zeros((nz,ny,nx));
    mask=np.tile(mask,(nz,1,1));
    #mask=np.transpose(mask, (2, 1, 0));

    n=0;
    if m2y==1:
         for year in xrange(fyear,lyear+1):
              for month in xrange(0,12):
                   n=n+months2days[month]
                   sdate="%4.4d%c%2.2d" % (year,datesep,month+1)
                   dnm=nc_read(prefix+sdate+'.nc','templvl');
                   templvl=templvl+np.squeeze(dnm.data)*mask*months2days[month]
                   dnm=nc_read(prefix+sdate+'.nc','salnlvl');
                   salnlvl=salnlvl+np.squeeze(dnm.data)*mask*months2days[month]
                   dnm=nc_read(prefix+sdate+'.nc','difisolvl');
                   difisolvl=difisolvl+10**(np.squeeze(dnm.data))*mask*months2days[month]
                   dnm=nc_read(prefix+sdate+'.nc','difdialvl');
                   difdialvl=difdialvl+10**(np.squeeze(dnm.data))*mask*months2days[month]
                   print sdate



def amoctime(args):
    ''' compute max amoc, amoc26 in time'''
    
    prefix,sdate = get_sdate_ini(root_folder) 
    # Latitude interval for maximum AMOC search
    lat=nc_read(prefix+sdate+'.nc','lat');
    lat1=20;
    lat2=60;
    lat3=26.5;
    ind1=np.min(np.where(lat>=lat1));
    ind2=np.max(np.where(lat<=lat2));
    ind3=np.max(np.where(lat<=lat3));
    amoc = {}

    for year in xrange(fyear,lyear+1):
        if m2y==1:
            for month in xrange(0,12):
#              n=n+months2days[month]
               sdate = "%4.4d%c%2.2d" % (year,datesep,month+1)
               dnm = nc_read(prefix+sdate+'.nc','mmflxd');
               dnm = np.copy(dnm[0,region-1,:,:])
               dnm[dnm>1e25] = np.nan
               amoc.setdefault('max',[]).append(np.nanmax(dnm[:,ind1-1:ind2-1])*1e-9)
               amoc.setdefault('26N',[]).append(np.nanmax(dnm[:,ind3-1])*1e-9)
               print sdate


        else:
            sdate = "%4.4d" % (year)
            dnm = nc_read(prefix+sdate+'.nc','mmflxd');
            dnm = np.copy(dnm[0,region-1,:,:])
            dnm[dnm>1e25] = np.nan
            amoc.setdefault('max',[]).append(np.nanmax(dnm[:,ind1-1:ind2-1])*1e-9)
            amoc.setdefault('26N',[]).append(np.nanmax(dnm[:,ind3])*1e-9)
            print sdate


    plt.plot(amoc['max'],'b',label='max')
    plt.plot(amoc['26N'],'r',label='26N')
    plt.legend(loc='lower right')
    return amoc


def amocmean(args):
    ''' compute mean amoc'''
    prefix,sdate = get_sdate_ini(root_folder) 
    # set initial values to zero
    amoc = set_ini_val_zero(prefix, sdate, 'x', 'y', 'depth')
    print 'mehmet', amoc.shape
    


    return amoc


def set_ini_val_zero(prefix,sdate,dim1,dim2,dim3):
    nx=ncgetdim(prefix+sdate+'.nc', dim1);
    ny=ncgetdim(prefix+sdate+'.nc', dim2);
    nz=ncgetdim(prefix+sdate+'.nc', dim3);
    var = np.zeros((nz,ny,nx))
    return var


# 1 for atlantic_arctic_ocean region
# 2 for indian_pacific_ocean region
# 3 for global_ocean
region = 1
amoc=amocmean(root_folder)
#amoc=amoctime(root_folder)
