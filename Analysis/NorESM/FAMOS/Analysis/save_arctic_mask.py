import numpy as np
import numpy.ma as ma
import sys
from netCDF4 import Dataset
#
sys.path.append('/Home/siv22/anu074/NorESM')
import NorESM_utils as utils
#read lon, lat info
filepath='/Data/skd/users/anu074/norstore/'
grid_1=Dataset(filepath+'grid_tnx1v2.nc')
lon_1=grid_1.variables['plon'][:]
lat_1=grid_1.variables['plat'][:]
#read data from one file to get the mask land-sea mask
y=100; m=1
data=Dataset(filepath+'NOIIA_T62_tn11_ctrl_default_02'+'/ocn/hist/'+'NOIIA_T62_tn11_ctrl_default_02'+'.micom.hm.'+str(y).zfill(4)+'-'+str(m).zfill(2)+'.nc')
sst=data.variables['sst'][:].squeeze()
mask=sst.mask
#pickup the ocean points
jinds,iinds=ma.where(1-mask)
#file where section indices are defined
secfile='/Data/skd/users/anu074/norstore/secindex_tnx1v1'
#then create the new Arctic mask
arctic_mask=utils.NorESM_masks('arctic', iinds, jinds, barents=True, secfile=secfile, grid='tripolar', lon=lon_1, lat=lat_1, baltic=False)
#
###############################
# SAVE THE NEW MASK TO NETCDF #
###############################
#
ylen,xlen=arctic_mask.shape
ncfile = Dataset('/work/anu074/NorESM_tnx1v2_arctic_mask.nc', 'w', format='NETCDF4')
#time = ncfile.createDimension('time', None)
y = ncfile.createDimension('y', ylen)
x = ncfile.createDimension('x', xlen)
#times = ncfile.createVariable('taxis','f8',('time',))
lat = ncfile.createVariable('lat','f4',('y','x'))
lon = ncfile.createVariable('lon','f4',('y','x'))
nc_var = ncfile.createVariable('mask','f4',('y','x',))
#put in the variable
lat[:]=lat_1[:]
lon[:]=lon_1[:]
nc_var[:]=arctic_mask[:]
#close the file
ncfile.close()