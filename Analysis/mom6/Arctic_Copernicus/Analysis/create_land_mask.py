import numpy as np


landmask_file = '/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/land_mask.nc'
gr1 = xr.open_dataset(landmask_file)

gr3 = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/ocean_hgrid.nc')
lon1 = gr3.x[1::2,1::2]
lat1 = gr3.y[1::2,1::2]

gr1['x'] = gr1.mask
gr1['y'] = gr1.mask
gr1

lon = np.copy(lon1)
dflon = xr.DataArray(lon, coords={'ny': gr1.ny, 'nx': gr1.nx},
             dims=['ny', 'nx'])
gr1['x'] = dflon

lat = np.copy(lat1)
dflat = xr.DataArray(lat, coords={'ny': gr1.ny, 'nx': gr1.nx},
             dims=['ny', 'nx'])
gr1['y'] = dflat

gr1.to_netcdf('/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/land_mask_v2.nc')
