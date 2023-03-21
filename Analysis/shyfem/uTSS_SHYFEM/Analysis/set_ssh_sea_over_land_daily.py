import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
from HCtFlood.kara import flood_kara
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
# import ESMF
from scipy.interpolate import interp1d
plt.ion()

root_folder = '/work/opa/mi19918/Projects/uTSS_SHYFEM/work/daily/'

# 'REG_40_41_0_u_velocity_0.01.nc'

# If false we will use original nemo grid 
variable = 'water_level'

# for tind in range(366,731): 
for tind in range(731,984): 
    print(tind)
    inputfile = root_folder + 'REG_' + str(tind) + '_' + str(tind+1) + '_0_' + variable + '_0.001.nc'
    ds = xr.open_dataset(inputfile)
    # extrapolate  ocean values into the land                         
    drowned_var = flood_kara(ds[variable], xdim='lon', ydim='lat')  
    ds1 = drowned_var.to_dataset(name=variable) 
    ds1 = ds1.isel(time=0)
    ds1 = ds1.drop('time')
    ds1 = ds1.isel(z=0)
    ds1 = ds1.drop('z')
    ds1['lon'] = ds['lon']  
    ds1['lat'] = ds['lat']  
    outputfile1 = root_folder + '/uTSS_OBC_daily_2021_' + variable + '_shyfem_levels_sea_over_land_high_res_' +  str(tind).zfill(3) + '.nc' 
    ds1.to_netcdf(outputfile1)




