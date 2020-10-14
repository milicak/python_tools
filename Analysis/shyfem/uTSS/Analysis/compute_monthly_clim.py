import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import xarray as xr

project_folder = '/work/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/'
root_folder = '/work/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/OUT/'

# list = sorted(glob.glob(root_folder+'*nos*'))
# df = xr.open_mfdataset(list, chunks={'time':5, 'node':20000})
# df = xr.open_mfdataset(list)

# salt = df['salinity'].groupby('time.month').mean('time')
# ds = salt.to_dataset(name='salinity')
# dnm = df.groupby('time.month').mean('time')
# outputfile = project_folder + 'uTSS_TS_monthly_clim.nc'
# dnm.to_netcdf(outputfile)


# list = sorted(glob.glob(root_folder+'*ous*'))
# df = xr.open_mfdataset(list, chunks={'time':5, 'node':20000})
# df = xr.open_mfdataset(list)
# dnm = df.groupby('time.month').mean('time')
# dnm.load()
# outputfile = project_folder + 'uTSS_UV_monthly_clim.nc'
# dnm.to_netcdf(outputfile)


list = sorted(glob.glob(root_folder+'*nos*'))
df = xr.open_mfdataset(list, chunks={'time':25, 'node':20000})
dnm = df.groupby('time.dayofyear').mean('time')
outputfile = project_folder + 'uTSS_TS_daily_clim.nc'
dnm.to_netcdf(outputfile)

# list = sorted(glob.glob(root_folder+'*ous*'))
# df = xr.open_mfdataset(list, chunks={'time':25, 'node':20000})
# dnm = df.groupby('time.dayofyear').mean('time')
# outputfile = project_folder + 'uTSS_UV_daily_clim.nc'
# dnm.to_netcdf(outputfile)
