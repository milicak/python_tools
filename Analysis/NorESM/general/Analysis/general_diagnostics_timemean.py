import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import sys
from netcdf_functions import nc_read
from netcdf_functions import ncgetdim

root_folder='/hexagon/work/milicak/archive/'

expid='NOIIA_T62_tn11_norems2_ctrl_tke'

fyear = 80; # first year
lyear = 100; # last year

m2y=1

# tripolar 1degree grid
grid_file='/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid.nc';
# tripolar 0.25degree grid
#grid_file = '/bcmhsm/milicak/RUNS/noresm/CORE2/Arctic/maps/grid_0_25degree.nc';
# bi-polar grid
#grid_file='/fimm/home/bjerknes/milicak/Analysis/NorESM/climatology/Analysis/grid_bipolar.nc';

mask=nc_read(grid_file,'pmask');

datesep='-';

if m2y==1:
    prefix=root_folder+expid+'/ocn/hist/'+expid+ '.micom.hm.'
else:
    prefix=root_folder+expid+'/ocn/hist/'+expid+ '.micom.hy.'


# Get dimensions and time attributes
if m2y==1:
    sdate="%4.4d%c%2.2d" % (fyear,datesep,1)
else:
    sdate="%4.4d" % (fyear)


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



