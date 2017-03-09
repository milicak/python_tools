import numpy as np
import sys
import os
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from netcdf_functions import nc_read
import matplotlib.pyplot as plt
#import netcdf
sys.path.insert(0,'/home/mil021/anaconda2/envs/esmpy/lib/python2.7/site-packages/')
import ESMF

plt.ion()


def create_locstream_spherical_16(coord_sys=ESMF.CoordSys.SPH_DEG, domask=False):
    """
    :param coord_sys: the coordinate system of the LocStream
    :param domask: a boolean to tell whether or not to add a mask
    :return: LocStream
    """
    if ESMF.pet_count() is not 1:
        raise ValueError("processor count must be 1 to use this function")


    lon_s4=np.array([17.6, 16.5, 16.05, 15.6, 15.1, 14.1, 13.0, 12.0, 10.0, 8.0, 4.0
                 , 4.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0
                 , 110.0, 120.0, 130.0, 140.0]);
    lat_s4=np.array([69.0, 70.6, 71.3, 72.02, 72.8, 73.8, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0
        , 81.0, 81.8, 81.8, 82.6, 83.0, 83.2, 83.1, 82.8, 82.5, 81.8, 79.7, 78.2
        , 78.2, 79.7]);

    locstream = ESMF.LocStream(lon_s4.shape[0], coord_sys=coord_sys)
    #locstream = ESMF.LocStream(16, coord_sys=coord_sys)
    deg_rad = np.pi
    if coord_sys == ESMF.CoordSys.SPH_DEG:
        deg_rad = 180

    locstram["ESMF:Lon"] = lon_s4
    locstram["ESMF:Lat"] = lat_s4
    #locstream["ESMF:Lon"] = [0.0, 0.5*deg_rad, 1.5*deg_rad, 2*deg_rad, 0.0, 0.5*deg_rad, 1.5*deg_rad, 2*deg_rad, 0.0, 0.5*deg_rad, 1.5*deg_rad, 2*deg_rad, 0.0, 0.5*deg_rad, 1.5*deg_rad, 2*deg_rad]
    #locstream["ESMF:Lat"] = [deg_rad/-2.0, deg_rad/-2.0, deg_rad/-2.0, deg_rad/-2.0, -0.25*deg_rad, -0.25*deg_rad, -0.25*deg_rad, -0.25*deg_rad, 0.25*deg_rad, 0.25*deg_rad, 0.25*deg_rad, 0.25*deg_rad, deg_rad/2.0, deg_rad/2.0, deg_rad/2.0, deg_rad/2.0]
    if domask:
        locstream["ESMF:Mask"] = np.array([1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=np.int32)

    return locstream




m = Basemap(width=8000000,height=8000000,
           resolution='l',projection='stere',
           lat_ts=40,lat_0=90,lon_0=0.)
m.drawcoastlines()
xpt,ypt = m(lon_s4,lat_s4)
m.plot(xpt,ypt,color='r')


grid1="/mnt/fimmhome/python_tools/Analysis/NorESM/APPLICATE/Analysis/noresm_ESMF_grid_tnx1v1.nc"
grid2="/mnt/fimmhome/Analysis/NorESM/climatology/Analysis/t00an1.nc"

# Create a uniform global latlon grid from a GRIDSPEC formatted file source grid
srcgrid = ESMF.Grid(filename=grid2,
                 filetype=ESMF.FileFormat.GRIDSPEC)

coord_sys=ESMF.CoordSys.SPH_DEG
domask=True


locstream = create_locstream_spherical_16(coord_sys=coord_sys, domask=domask)

srcfield = ESMF.Field(srcgrid, staggerloc=ESMF.StaggerLoc.CENTER)
tempwoa = nc_read(grid2,'t')
tempwoa = tempwoa[0,:,:]
srcfield.data[:] = np.transpose(tempwoa)

# create a field on the locstream
dstfield = ESMF.Field(locstream, name='dstfield')

# create an object to regrid data from the source to the destination field
dst_mask_values=None
if domask:
        dst_mask_values=np.array([0])

regrid = ESMF.Regrid(srcfield, dstfield,
                    regrid_method=ESMF.RegridMethod.BILINEAR,
                    unmapped_action=ESMF.UnmappedAction.IGNORE,
                    dst_mask_values=dst_mask_values)

# do the regridding from source to destination field
dstfield = regrid(srcfield, dstfield)


sys.exit()






# Create a field on the centers of the grid
#field1 = ESMF.Field(srcgrid, staggerloc=ESMF.StaggerLoc.CENTER,ndbounds=[33, 1])    # level=33 time =1


# read the variable
fname = netcdf.netcdf_file(grid2, 'r')
salt = fname.variables['s']
salt=salt[0,:,:]
field1.data[:] = np.transpose(salt)


# Read a field from the file
field1.read(filename=grid2, variable="s")

# however you want the salt to be read is up to you
field1.data[:] = salt

# Create a uniform global latlon grid from a SCRIP formatted file
dstgrid = ESMF.Grid(filename=grid1, filetype=ESMF.FileFormat.SCRIP)

#dstgrid = ESMF.Grid(filename=grid1, filetype=ESMF.FileFormat.SCRIP,
#                 add_corner_stagger=True)



# Create a field on the centers of the grid
field2 = ESMF.Field(dstgrid, staggerloc=ESMF.StaggerLoc.CENTER)


field2.data[:] = 0


# Set up a regridding object between source and destination
regridS2D = ESMF.Regrid(field1, field2,
                        regrid_method=ESMF.RegridMethod.BILINEAR)

field2 = ESMF.Regrid(field1, field2)

