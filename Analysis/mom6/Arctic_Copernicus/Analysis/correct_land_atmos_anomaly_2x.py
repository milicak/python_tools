import numpy as np


df = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/2X/atmos_mosaic_tile1Xocean_mosaic_tile1.nc')
ind=df.xgrid_area.argmax()
df.xgrid_area[ind] = df.xgrid_area[ind+1]
df.to_netcdf('/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/2X/atmos_mosaic_tile1Xocean_mosaic_tile1iv2.nc')

df=xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/2X/land_mosaic_tile1Xocean_mosaic_tile1.nc')
ind=df.xgrid_area.argmax()
df.xgrid_area[ind] = df.xgrid_area[ind+1]
df.to_netcdf('/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/2X/land_mosaic_tile1Xocean_mosaic_tile1v2.nc')
df.close()
del df


df = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus/2X/atmos_mosaic_tile1Xland_mosaic_tile1.nc')


