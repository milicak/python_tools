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


#m = Basemap(width=8000000,height=8000000,
#           resolution='l',projection='stere',
#           lat_ts=40,lat_0=90,lon_0=0.)
#m.drawcoastlines()


#grid1 = "/mnt/fimmhome/python_tools/Analysis/NorESM/APPLICATE/Analysis/noresm_ESMF_grid_tnx1v1.nc"
grid1 = "/mnt/fimmhome/python_tools/Analysis/NorESM/APPLICATE/Analysis/noresm_ESMF_grid_tnx1v1_nohalo.nc"
grid2 = "/mnt/fimmhome/Analysis/NorESM/climatology/Analysis/t00an1.nc"
grd="/mnt/fimmhome/Analysis/NorESM/climatology/Analysis/grid.nc"

##### From NorESM to WOA Grid ###############
# Create a uniform global latlon grid from a GRIDSPEC formatted file source grid
srcgrid = ESMF.Grid(filename=grid1,
                 filetype=ESMF.FileFormat.SCRIP)

# Create a field on the centers of the grid
field1 = ESMF.Field(srcgrid, staggerloc=ESMF.StaggerLoc.CENTER)
datafile = "/mnt/hexwork/archive/NBF1850_f19_tn11_sst_sss_rlx_01/ocn/hist/NBF1850_f19_tn11_sst_sss_rlx_01.micom.hy.0070.nc"

tempsrc = nc_read(datafile,'templvl')
tempsrc = tempsrc[0,0,:-1,:]
field1.data[:] = np.transpose(tempsrc)

dstgrid = ESMF.Grid(filename=grid2, filetype=ESMF.FileFormat.GRIDSPEC)
# Create a field on the centers of the grid
field2 = ESMF.Field(dstgrid, staggerloc=ESMF.StaggerLoc.CENTER)

field2.data[:] = 0

# Set up a regridding object between source and destination
regridS2D = ESMF.Regrid(field1, field2,
                        regrid_method=ESMF.RegridMethod.BILINEAR)
field2 = regridS2D(field1, field2)

sys.exit()


##### From NorESM to WOA Grid ###############

srcgrid = ESMF.Grid(filename=grd,filetype=ESMF.FileFormat.GRIDSPEC,coord_names=["vlon",
                                                                      "vlat"])
field1 = ESMF.Field(srcgrid, staggerloc=ESMF.StaggerLoc.CENTER)
field1.data[:]=np.transpose(tempsrc)
dstgrid = ESMF.Grid(filename=grid2, filetype=ESMF.FileFormat.GRIDSPEC)
field2 = ESMF.Field(dstgrid, staggerloc=ESMF.StaggerLoc.CENTER)
field2.data[:] = 0
regridS2D = ESMF.Regrid(field1,
                        field2,regrid_method=ESMF.RegridMethod.BILINEAR)
field2 = regridS2D(field1, field2)

sys.exit()

##### From WOA to NorESM Grid ###############

# Create a uniform global latlon grid from a GRIDSPEC formatted file source grid
srcgrid = ESMF.Grid(filename=grid2,
                 filetype=ESMF.FileFormat.GRIDSPEC)

# Create a field on the centers of the grid
field1 = ESMF.Field(srcgrid, staggerloc=ESMF.StaggerLoc.CENTER)

# read the variable
#fname = netcdf.netcdf_file(grid2, 'r')
#salt = fname.variables['t']
#salt=salt[0,:,:]
#field1.data[:] = np.transpose(salt)

# Read a field from the file 2nd way
#field1.read(filename=grid2, variable="s")

# however you want the salt to be read is up to you
tempwoa = nc_read(grid2,'t')
field1.data[:] = np.transpose(tempwoa[0,:,:])

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
field2 = regridS2D(field1, field2)

datafile = "/mnt/hexwork/archive/NBF1850_f19_tn11_sst_sss_rlx_01/ocn/hist/NBF1850_f19_tn11_sst_sss_rlx_01.micom.hy.0070.nc"

tempsrc = nc_read(datafile,'templvl')
tempsrc = tempsrc[0,0,:,:]
sys.exit()

##############################################




field1 = ESMF.Field(srcgrid, staggerloc=ESMF.StaggerLoc.CENTER)


sys.exit()


srcgrid = ESMF.Grid(filename=grd,filetype=ESMF.FileFormat.GRIDSPEC,coord_names=["plon",
                                                                      "plat"])
field1 = ESMF.Field(srcgrid, staggerloc=ESMF.StaggerLoc.CENTER)
field1.data[:]=np.transpose(tempsrc)
dstgrid = ESMF.Grid(filename=grid2, filetype=ESMF.FileFormat.GRIDSPEC)
field2 = ESMF.Field(dstgrid, staggerloc=ESMF.StaggerLoc.CENTER)
field2.data[:] = 0
