import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
from HCtFlood.kara import flood_kara
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import ESMF
from mpl_toolkits.basemap import Basemap                                            
from scipy.interpolate import interp1d
plt.ion()

root_folder  = '/work/opa/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/monthly_clim/'


# If false we will use original nemo grid 
variable = 'v_velocity'
# variable = 'u_velocity'

for tind in range(0,12): 
    print(tind)
    inputfile = root_folder + '/uTSS_OBC_monthly_' + variable + '_nemo_levels_low_res_' +  np.str(tind).zfill(2) + '.nc' 
    # inputfile = root_folder + '/uTSS_OBC_monthly_' + variable + '_nemo_levels_' +  np.str(tind).zfill(2) + '.nc' 
    ds = xr.open_dataset(inputfile)
    # extrapolate salinity ocean values into the land                         
    drowned_var = flood_kara(ds[variable], xdim='x', ydim='y', zdim='z')  
    ds1 = drowned_var.to_dataset(name=variable) 
    ds1 = ds1.isel(time=0)
    ds1 = ds1.drop('time')
    ds1['lon'] = ds['lon']  
    ds1['lat'] = ds['lat']  
    # extrapolate vertically to fill NaN values                               
    ds1 = ds1.ffill(dim='z')                           
    # outputfile1 = root_folder + '/uTSS_OBC_monthly_' + variable + '_nemo_levels_sea_over_land_' +  np.str(tind).zfill(2) + '.nc' 
    outputfile1 = root_folder + '/uTSS_OBC_monthly_' + variable + '_nemo_levels_sea_over_land_low_res_' +  np.str(tind).zfill(2) + '.nc' 
    ds1.to_netcdf(outputfile1)




