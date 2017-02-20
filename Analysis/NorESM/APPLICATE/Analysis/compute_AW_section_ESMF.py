import numpy as np
import sys
import os
#from netcdf_functions import nc_read
import matplotlib.pyplot as plt
sys.path.insert(0,'/home/mil021/anaconda2/envs/esmpy/lib/python2.7/site-packages/')
import ESMF
#import netcdf

plt.ion()


grid1="/mnt/fimmhome/python_tools/Analysis/NorESM/APPLICATE/Analysis/noresm_ESMF_grid_tnx1v1.nc"
grid2="/mnt/fimmhome/Analysis/NorESM/climatology/Analysis/s00an1.nc"

# Create a uniform global latlon grid from a GRIDSPEC formatted file source grid
srcgrid = ESMF.Grid(filename=grid2,
                 filetype=ESMF.FileFormat.GRIDSPEC)
sys.exit()

# Create a field on the centers of the grid
#field1 = ESMF.Field(srcgrid, staggerloc=ESMF.StaggerLoc.CENTER,ndbounds=[33, 1])    # level=33 time =1
field1 = ESMF.Field(srcgrid, staggerloc=ESMF.StaggerLoc.CENTER)


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

