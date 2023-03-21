import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import glob

gr = xr.open_dataset('/data/opa/bs-mod/upstream_bs-nrt/geodta/mesh_mask_bs-nrt_s5_smth_changeBosp_Marm_sill_nemo.nc')

root_folder = '/work/opa/mi19918/Projects/uTSS_SHYFEM/work/daily/'
znemo = np.transpose(np.copy(gr.gdept_1d))   
znemo = znemo[:,0]

# ls1 = sorted(glob.glob(root_folder+'uTSS_OBC_daily_2021_u_velocity_shyfem_levels_sea_over_land_low_res*'))
# ls1 = ls1[731:]
# for fname in ls1:
#     print(fname)
#     ds = xr.open_dataset(fname)
#     dsnemo  = ds.interp(level=znemo)    
#     dsnemo = dsnemo.bfill('level')
#     outputfile = fname[0:81] + 'nemo_levels_' + fname[-6:-3] + '.nc'
#     dsnemo.to_netcdf(outputfile)
#     dsnemo.close()
#
#
# ls1 = sorted(glob.glob(root_folder+'uTSS_OBC_daily_2021_v_velocity_shyfem_levels_sea_over_land_low_res*'))
# ls1 = ls1[731:]
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
ls1 = ls1[731:]
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
ls1 = ls1[731:]
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
