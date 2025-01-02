import numpy as np

df = xr.open_dataset('~/dataset/MOM6/Arctic_Copernicus/ocean_topog.nc')

ds = df.depth[::2,::2]
ds = ds.to_dataset(name='depth')

ds.to_netcdf('/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/2X/ocean_topog.nc')
ds.close()
