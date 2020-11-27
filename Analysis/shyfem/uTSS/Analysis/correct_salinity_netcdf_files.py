import numpy as np
import os
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
from HCtFlood.kara import flood_kara
import matplotlib.pyplot as plt

root_folder = '/archive/milicak/shyfem/uTSS/Exp_2016_analysis_newTSIC/OUT/'
top_coef = (1-0.09*(1-np.exp(-(30-23)/(20))**10))


# for tind in range(0,12):
for tind in range(0,1460):
    print(tind)
    filename = root_folder + '/uTSS_lobc_chunk_' +  np.str(tind).zfill(4) + '.nos.nc'
    df = xr.open_dataset(filename)
    dnm = df.salinity[0,:,:]
    meets_cond3 = (dnm>30)
    meets_cond4 = (dnm<20)
    meets_cond2 = (dnm>=28) & (dnm<=30)
    meets_cond = (dnm>23) & (dnm<28)
    dnm.values[meets_cond3]=dnm.values[meets_cond3]*0.999998
    dnm.values[meets_cond4]=dnm.values[meets_cond4]*1.000001
    dnm.values[meets_cond2]=dnm.values[meets_cond2]*top_coef
    dnm.values[meets_cond]=dnm.values[meets_cond]*(1-0.09*(1-np.exp(-(dnm.values[meets_cond]-23)/(20))**10))
    df.salinity[0,:,:] = dnm
    df.to_netcdf('dnm.nc')
    cmnd = 'mv dnm.nc ' + filename
    os.system(cmnd)


