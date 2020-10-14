import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import xarray as xr

project_folder = '/work/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/'
root_folder = '/work/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/OUT/'

# for ind in range(1,366):
#     print(ind)
#     files = ind + np.array([365*1, 365*2, 365*3]) 
#     list = sorted(glob.glob(root_folder+ '*' + np.str(files[0]).zfill(4) + '*nos*'))
#     list.extend(sorted(glob.glob(root_folder+ '*' + np.str(files[1]).zfill(4) + '*nos*')))
#     list.extend(sorted(glob.glob(root_folder+ '*' + np.str(files[2]).zfill(4) + '*nos*')))
#     df = xr.open_mfdataset(list)
#     ds = df.mean('time')
#     # dnm = df.groupby('time.dayofyear').mean('time')
#     outputfile = project_folder + 'daily_clim/'  + 'uTSS_TS_daily_clim_' + np.str(ind).zfill(4) + '.nc'
#     ds.to_netcdf(outputfile)



for ind in range(1,366):
    print(ind)
    files = ind + np.array([365*1, 365*2, 365*3]) 
    list = sorted(glob.glob(root_folder+ '*' + np.str(files[0]).zfill(4) + '*ous*'))
    list.extend(sorted(glob.glob(root_folder+ '*' + np.str(files[1]).zfill(4) + '*ous*')))
    list.extend(sorted(glob.glob(root_folder+ '*' + np.str(files[2]).zfill(4) + '*ous*')))
    df = xr.open_mfdataset(list)
    ds = df.mean('time')
    # dnm = df.groupby('time.dayofyear').mean('time')
    outputfile = project_folder + 'daily_clim/'  + 'uTSS_UV_daily_clim_' + np.str(ind).zfill(4) + '.nc'
    ds.to_netcdf(outputfile)

