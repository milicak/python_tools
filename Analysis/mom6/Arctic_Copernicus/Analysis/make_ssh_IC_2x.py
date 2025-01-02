import numpy as np

df = xr.open_dataset('~/dataset/MOM6/Arctic_Copernicus/Glorys_ssh_IC.nc')

ds = df.ssh[:,::2,::2]
ds = ds.to_dataset(name='ssh')
ds['lon'] = df.lon[::2,::2]
ds['lat'] = df.lat[::2,::2]

ds.to_netcdf('/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/2X/Glorys_ssh_IC.nc',unlimited_dims=['time'])
ds.close()
