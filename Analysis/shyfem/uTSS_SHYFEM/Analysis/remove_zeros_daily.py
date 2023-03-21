import numpy as np
import glob


ls1 = sorted(glob.glob('/work/opa/mi19918/Projects/uTSS_SHYFEM/work/OUT/*nos*'))
# ls1 = ls1[791:]
for fname in ls1:
    print(fname)
    df = xr.open_dataset(fname)
    df = df.where(df!=0)
    df.to_netcdf(fname)
    df.close; del df



ls1 = sorted(glob.glob('/work/opa/mi19918/Projects/uTSS_SHYFEM/work/OUT/*ous*'))
for fname in ls1:
    print(fname)
    df = xr.open_dataset(fname)
    df = df.where(df!=0)
    df.to_netcdf(fname)
    df.close; del df
