import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import sys
def nc_read(fname,variable):
    # how to use this subroutine is from netcdf_functions import nc_read
    ncfile = Dataset(fname, 'r', format='NETCDF4')
    list(ncfile.variables)
    v1=ncfile.variables[variable][:].copy();
    return v1

def ncgetdim(fname,variable):
    # how to use this subroutine is from netcdf_functions import nc_read
    ncfile = Dataset(fname, 'r', format='NETCDF4')
    v1=len(ncfile.dimensions[variable]);
    return v1
