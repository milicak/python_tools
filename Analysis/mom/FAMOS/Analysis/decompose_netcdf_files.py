# first goto /export directory run there
import numpy as np
import netCDF4
import os
fname = '00010101.ocean_month.nc'

dset = netCDF4.Dataset(fname)
varnames = dset.variables.keys()

for varname in varnames:
    print varname
    outname = varname+'_'+fname
    netcdfcommand = "ncks -v "+varname+" "+fname+" "+outname
    os.system(netcdfcommand)

