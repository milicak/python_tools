import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import glob
import os

time = pd.date_range("2020-01-01", freq="D", periods=365 +366)
root_folder = '/work/opa/mi19918/Projects/uTSS_SHYFEM/work/daily/'
ls1 = sorted(glob.glob(root_folder+'uTSS_OBC_daily_2021_salinity_nemo_levels*'))
for ind,fname in enumerate(ls1):
    print(fname)
    fnamenew = fname[:-6] + str(time[ind])[:10] + '.nc'
    cmnd = 'mv ' + fname + '  ' + fnamenew
    os.system(cmnd)

# ls1 = sorted(glob.glob(root_folder+'uTSS_OBC_daily_2021_u_velocity_shyfem_levels_sea_over_land_low_res*'))
# ls1 = ls1[0:366]
# for fname in ls1:
#     print(fname)
#     ds = xr.open_dataset(fname)
#     dsnemo  = ds.interp(level=znemo)    
#     dsnemo = dsnemo.bfill('level')
#     outputfile = fname[0:81] + 'nemo_levels_' + fname[-6:-3] + '.nc'
#     dsnemo.to_netcdf(outputfile)
#     dsnemo.close()
#

# ls1 = sorted(glob.glob(root_folder+'uTSS_OBC_daily_2021_v_velocity_shyfem_levels_sea_over_land_low_res*'))
# ls1 = ls1[0:366]
# for fname in ls1:
#     print(fname)
#     ds = xr.open_dataset(fname)
#     dsnemo  = ds.interp(level=znemo)    
#     dsnemo = dsnemo.bfill('level')
#     outputfile = fname[0:81] + 'nemo_levels_' + fname[-6:-3] + '.nc'
#     dsnemo.to_netcdf(outputfile)
#     dsnemo.close()
#

ls1 = sorted(glob.glob(root_folder+'uTSS_OBC_daily_2021_temperature_shyfem_levels_sea_over_land_high_res*'))
ls1 = ls1[0:366]
for fname in ls1:
    print(fname)
    ds = xr.open_dataset(fname)
    dsnemo  = ds.interp(level=znemo)    
    dsnemo = dsnemo.bfill('level')
    dsnemo = dsnemo.ffill('level')
    outputfile = fname[0:81] + '_nemo_levels_' + fname[-6:-3] + '.nc'
    dsnemo.to_netcdf(outputfile)
    dsnemo.close()

ls1 = sorted(glob.glob(root_folder+'uTSS_OBC_daily_2021_salinity_shyfem_levels_sea_over_land_high_res*'))
ls1 = ls1[0:366]
for fname in ls1:
    print(fname)
    ds = xr.open_dataset(fname)
    dsnemo  = ds.interp(level=znemo)    
    dsnemo = dsnemo.bfill('level')
    dsnemo = dsnemo.ffill('level')
    outputfile = fname[0:78] + '_nemo_levels_' + fname[-6:-3] + '.nc'
    dsnemo.to_netcdf(outputfile)
    dsnemo.close()

# old code
# variable = 'temperature'
# for tind in range(366,731): 
#     print(tind)
#     inputfile = root_folder + 'REG_' + str(tind) + '_' + str(tind+1) + '_0_' + variable + '_0.001.nc'
#     ds = xr.open_dataset(inputfile)
#     dsnemo  = ds.interp(level=znemo)    
#     dsnemo = dsnemo.bfill('level')
#     outputfile = fname[0:70] + 'temperature_nemo_levels_' + str(tind).zfill(3) + '.nc'
#     dsnemo.to_netcdf(outputfile)
#     dsnemo.close()
#
#
# variable = 'salinity'
# for tind in range(366,731): 
#     print(tind)
#     inputfile = root_folder + 'REG_' + str(tind) + '_' + str(tind+1) + '_0_' + variable + '_0.001.nc'
#     ds = xr.open_dataset(inputfile)
#     dsnemo  = ds.interp(level=znemo)    
#     dsnemo = dsnemo.bfill('level')
#     outputfile = fname[0:70] + 'salinity_nemo_levels_' + str(tind).zfill(3) + '.nc'
#     dsnemo.to_netcdf(outputfile)
#     dsnemo.close()
#
