import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd

rho_0 = 1025
Cp = 3996

# root_folder = '/archive2/milicak/mitgcm/sose/'
root_folder = '/archive/milicak/MITgcm/Projects/sose/'
expid = 'Exp01_0'
print(expid)
outname = root_folder + expid + '/' + expid + '_meridional_heat_transport.nc'

fnames = root_folder + expid + '/*ADVy_TH_*'

list = sorted(glob.glob(fnames))

df = xr.open_mfdataset(list)
df = df.sel(time=slice('2007','2012'))
df = df.mean('time')
df = df.sum(('i','k'))
ds = rho_0*Cp*df

ds.to_netcdf(outname)


