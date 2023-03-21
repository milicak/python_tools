import numpy as np


ds=df.temperature[0,:,:]
ds=ds.where(ds!=0)
dnm=ds.notnull().sum('level')
ds.isel(level=dnm-1)
