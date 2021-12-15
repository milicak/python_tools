import numpy as np
import xarray as xr

df = xr.open_dataset('/work/opa/mi19918/Projects/mom6/Arctic_Copernicus/INPUT/NEMO_TS_IC.nc')

df['temp'] = df.temp.where(df.temp>-1000)
df['salt'] = df.salt.where(df.salt>-1000)
df['temp'] = df.temp.ffill(dim='zt') 
df['salt'] = df.salt.ffill(dim='zt') 

df.to_netcdf('/work/opa/mi19918/Projects/mom6/Arctic_Copernicus/INPUT/NEMO_TS_IC.nc')
