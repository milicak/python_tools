import numpy as np
import sys
import os
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib.pyplot as plt
#import netcdf
#sys.path.insert(0,'/home/mil021/anaconda2/envs/esmpy/lib/python2.7/site-packages/')
import ESMF

plt.ion()


#m = Basemap(width=8000000,height=8000000,
#           resolution='l',projection='stere',
#           lat_ts=40,lat_0=90,lon_0=0.)
#m.drawcoastlines()

# 0.5 degree
# grid1 = 'mom6_05deg_esmf_meshinfo.nc'
grid1 = 'mom6_0666deg_esmf_meshinfo.nc'
grid2 = 'zos_AVISO_L4_199210-201012.nc'

##### From WOA to NorESM Grid ###############

# create source grid from WOA grid/data netcdf file
srcgrid = ESMF.Grid(filename=grid2, filetype=ESMF.FileFormat.GRIDSPEC)

# Create a uniform global latlon grid from a GRIDSPEC formatted file source grid
# create a destination grid file
dstgrid = ESMF.Grid(filename=grid1,filetype=ESMF.FileFormat.GRIDSPEC,
                    coord_names=["plon","plat"])
# if variables of lon and lat are lon and lat
#distgrid = ESMF.Grid(filename=grid1,filetype=ESMF.FileFormat.SCRIP)

# Create a field on the centers of the grid
field1 = ESMF.Field(srcgrid, staggerloc=ESMF.StaggerLoc.CENTER)
tempsrc = xr.open_dataset(grid2,decode_times=False)
woadata = []

# create a temporary source file for the surface size of tempsrc = [1,102,720,1440]
tmpsrc = np.copy(tempsrc['zos'].data[3:-12,:,:])
tmpsrc[tmpsrc>10] = np.nan
tmpsrc = np.nanmean(tmpsrc, axis=0)

# set it to zero
field1.data[:] = 0
field1.data[:] = np.transpose(tmpsrc)

# Create a field on the centers of the grid
field2 = ESMF.Field(dstgrid, staggerloc=ESMF.StaggerLoc.CENTER)
# set it to zero
field2.data[:] = 0

# Set up a regridding object between source and destination
regridS2D = ESMF.Regrid(field1, field2,
                       regrid_method=ESMF.RegridMethod.BILINEAR)
field2 = regridS2D(field1, field2)
dnm = np.expand_dims(field2.data, axis=2)
woadata = dnm[:,:,0]


# create data fram for the field2
#df = xr.DataArray(data=field2.data,dims=['x','y'],name='sstwoa_noresm')
#df = xr.DataArray(data=woadata,dims=['x','y','depth'],name='tempwoa_noresm')
df = xr.DataArray(data=woadata,dims=['x','y'],name='ssh_Aviso_mom6')

# df.to_netcdf('SSH_Aviso_mom6_0_5deg.nc')
df.to_netcdf('SSH_Aviso_mom6_0_666deg.nc')

