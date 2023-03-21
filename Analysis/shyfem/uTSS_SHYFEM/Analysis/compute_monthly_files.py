import numpy as np
import glob

root_folder = '/work/opa/mi19918/Projects/uTSS_SHYFEM/work/OUT/'

ls1 = sorted(glob.glob(root_folder+'*nos*'))
# ls1 = sorted(glob.glob(root_folder+'*ous*'))

# first year 2020
# ls1 = ls1[0:366]
# second year 2021
ls1 = ls1[366:366+365]
# for 2020
df1 = xr.open_mfdataset(ls1[1:])
# for 2021
df1 = xr.open_mfdataset(ls1)
df = df1.groupby('time.month').mean('time')

df.to_netcdf('/work/opa/mi19918/Projects/uTSS_SHYFEM/work/monthly/uTSS_lobc_chunk_monthly_2021.nos.nc')
# df.to_netcdf('/work/opa/mi19918/Projects/uTSS_SHYFEM/work/monthly/uTSS_lobc_chunk_monthly_2021.ous.nc')
