import xmitgcm
import numpy as np
import numpy.ma as ma
import glob
import sys
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LinearSegmentedColormap
import os
import xarray as xr


# compression options
comp = dict(zlib=True, complevel=5)

root_folder = '/work/opa/mi19918/Projects/nemo/BS/'
expid = 'BS-NRT_MI_2.3_fcorv2'

folderpath = root_folder + expid + '/rebuilt/'
cmd1 = folderpath + '*.nc'

ls1 = sorted(glob.glob(cmd1))

for ind in range(0,len(ls1)):
    fname = ls1[ind]
    print(fname)
    df = xr.open_dataset(fname)
    outname = folderpath+'dnm.nc'
    if not os.path.isfile(outname):
        encoding = {var: comp for var in df.data_vars}
        if bool(encoding):
            df.to_netcdf(outname, encoding=encoding)
            cmd2 = 'mv ' + outname + ' ' + fname
            os.system(cmd2)



