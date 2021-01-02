import os
import numpy as np
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import matplotlib.pyplot as plt
os.environ['NUMPY_EXPERIMENTAL_ARRAY_FUNCTION'] = '0'

from matplotlib import pyplot as plt
import numpy as np
import xarray as xr
import gsw
from xhistogram.xarray import histogram

root_folder = '/archive2/milicak/mitgcm/sose/'
project_name = 'Exp03_0'
fname = root_folder + project_name + '/grid.nc'
gr  = xr.open_dataset(fname)

list = sorted(glob.glob('/archive2/milicak/mitgcm/sose/Exp01_0/THETA*.nc'))
comp = dict(zlib=True, complevel=5)

# for ind in np.arange(73,len(list)):
for ind in np.arange(0,len(list)):
    inname = list[ind][-13:]
    fname = root_folder + project_name + '/' + 'THETA_' + inname
    dft = xr.open_dataset(fname)
    fname = root_folder + project_name + '/' + 'SALT_' + inname
    dfs = xr.open_dataset(fname)

    dft = dft.rename_dims({'i': 'XC', 'j': 'YC', 'i_g': 'XG', 'j_g': 'YG', 'k': 'Z', 'k_l': 'Zl', 'k_p1': 'Zp1', 'k_u': 'Zu'})
    dfs = dfs.rename_dims({'i': 'XC', 'j': 'YC', 'i_g': 'XG', 'j_g': 'YG', 'k': 'Z', 'k_l': 'Zl', 'k_p1': 'Zp1', 'k_u': 'Zu'})

    dfs1 = dfs.isel(time=0)
    dft1 = dft.isel(time=0)
    sigma2 = xr.apply_ufunc(gsw.sigma2, dfs1.SALT, dft1.THETA,
                                dask='parallelized', output_dtypes=[dfs1.SALT.dtype])
    df1 = sigma2.to_dataset(name='sigma2')
    fname = root_folder + project_name + '/' + 'SIGMA2_' + inname
    encoding = {var: comp for var in df1.data_vars}
    print(fname)
    df1.to_netcdf(fname, encoding=encoding)




