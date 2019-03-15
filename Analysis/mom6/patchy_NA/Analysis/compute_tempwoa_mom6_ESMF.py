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
grid1 = 'mom6_05deg_esmf_meshinfo.nc'
datagrid = xr.open_dataset('analysis_vgrid_lev35.v1.nc')

grid2 = '/cluster/NS2345K/noresm_diagnostics/packages/MICOM_DIAG/obs_data/WOA13/0.25deg/woa13_decav_t00_04.nc'

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
for depthind in range(0,102):
    print(depthind)
    tmpsrc = np.copy(tempsrc['t_an'].data[0,depthind,:,:])
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
    if depthind == 0:
        woadata = dnm
    else:
        woadata = np.append(woadata,dnm,axis=2);



nx = woadata.shape[0]
ny = woadata.shape[1]
dnm = np.zeros([nx,ny,datagrid.zt.shape[0]])
for ii in range(0,nx):
    for jj in range(0,ny):
        dnm[ii,jj,:] = np.interp(datagrid.zt,tempsrc.depth,woadata[ii,jj,:])


# create data fram for the field2
#df = xr.DataArray(data=field2.data,dims=['x','y'],name='sstwoa_noresm')
#df = xr.DataArray(data=woadata,dims=['x','y','depth'],name='tempwoa_noresm')
df = xr.DataArray(data=dnm,dims=['x','y','depth'],name='tempwoa_orca')

df.to_netcdf('tempwoa_mom6_0_5deg.nc')

