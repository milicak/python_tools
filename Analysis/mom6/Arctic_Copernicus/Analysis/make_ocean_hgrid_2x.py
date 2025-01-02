import numpy as np

df = xr.open_dataset('~/dataset/MOM6/Arctic_Copernicus/ocean_hgrid_lon-100_260.nc')

ds = df.x[::2,::2]
ds = ds.to_dataset(name='x')
ds['y'] = df.y[::2,::2]
ds['angle_dx'] = df.angle_dx[::2,::2]*180.0/np.pi
ds['dx'] = df.dx[::2,::2]*2
ds['dy'] = df.dy[::2,::2]*2
ds['area'] = df.area[::2,::2]*4
ds['tile'] = df.tile

ds.to_netcdf('/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/2X/ocean_hgrid.nc')
ds.close()
