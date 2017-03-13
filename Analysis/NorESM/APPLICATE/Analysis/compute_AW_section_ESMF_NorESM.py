import numpy as np
import sys
import os
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from netcdf_functions import nc_read
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
#import netcdf
sys.path.insert(0,'/home/mil021/anaconda2/envs/esmpy/lib/python2.7/site-packages/')
import ESMF

plt.ion()

# section locations
lon_s4=np.array([17.6, 16.5, 16.05, 15.6, 15.1, 14.1, 13.0, 12.0, 10.0, 8.0, 4.0
             , 4.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0
             , 110.0, 120.0, 130.0, 140.0]);
lat_s4=np.array([69.0, 70.6, 71.3, 72.02, 72.8, 73.8, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0
    , 81.0, 81.8, 81.8, 82.6, 83.0, 83.2, 83.1, 82.8, 82.5, 81.8, 79.7, 78.2
    , 78.7, 79.7]);

x = np.linspace(1, lon_s4.shape[0], num=lon_s4.shape[0], endpoint=True)
f = interp1d(x,lon_s4)
g = interp1d(x,lat_s4)
xnew = np.linspace(1, lon_s4.shape[0], num=2*lon_s4.shape[0],
                   endpoint=True)

lon_s4new = f(xnew)
lat_s4new = g(xnew)

coord_sys=ESMF.CoordSys.SPH_DEG
domask=True
# create locstream
locstream = ESMF.LocStream(lon_s4new.shape[0], name="Atlantic Inflow Section", coord_sys=coord_sys)
# appoint the section locations
locstream["ESMF:Lon"] = lon_s4new
locstream["ESMF:Lat"] = lat_s4new
if domask:
    locstream["ESMF:Mask"] = np.array(np.ones(lon_s4new.shape[0]), dtype=np.int32)

plt.figure()
m = Basemap(width=8000000,height=8000000,
           resolution='l',projection='stere',
           lat_ts=40,lat_0=90,lon_0=0.)
m.drawcoastlines()
xpt,ypt = m(lon_s4new,lat_s4new)
m.plot(xpt,ypt,'-o',color='b')

# The data file should be in global latlon grid from a GRIDSPEC formatted file source grid
datafname = "/mnt/fimmexport/bckup/noresm/CORE2/Arctic/DATA/NorESM/NOIIA_T62_tn11_sr10m60d_01_temperature_pendatal_1-300.nc"
gridfile = "/mnt/fimmhome/python_tools/Analysis/NorESM/APPLICATE/Analysis/noresm_ESMF_grid_tnx1v1_nohalo.nc"
tempwoa = nc_read(datafname,'temp')
zt = nc_read(datafname,'depth')
lon=nc_read(datafname,'tlon')
lat=nc_read(datafname,'tlat')

sys.exit()

secfield = np.zeros((lon_s4new.shape[0],tempwoa.shape[0]))

for kind in range(0,tempwoa.shape[0]):
    print 'indice = ', kind
    # Create a uniform global latlon grid from a GRIDSPEC formatted file source grid
    srcgrid = ESMF.Grid(filename=gridfile,
                     filetype=ESMF.FileFormat.GRIDSPEC)

    srcfield = ESMF.Field(srcgrid, staggerloc=ESMF.StaggerLoc.CENTER)
    temp = tempwoa[kind,:,:]
    srcfield.data[:] = np.transpose(temp)

    # create a field on the locstream
    dstfield = ESMF.Field(locstream, name='dstfield')
    dstfield.data[:] = 0.0

    # create an object to regrid data from the source to the destination field
    dst_mask_values=None
    if domask:
            dst_mask_values=np.array([0])

    regrid = ESMF.Regrid(srcfield, dstfield,
                        #regrid_method=ESMF.RegridMethod.NEAREST_STOD,
                        regrid_method=ESMF.RegridMethod.BILINEAR,
                        unmapped_action=ESMF.UnmappedAction.IGNORE,
                        dst_mask_values=dst_mask_values)

    # do the regridding from source to destination field
    dstfield = regrid(srcfield, dstfield)
    secfield[:,kind] = dstfield.data


plt.figure()
plt.pcolor(xnew,-zt,np.transpose(secfield),vmin=-1.5,vmax=7);plt.colorbar()



