import numpy as np
import xesmf as xe

WOA = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/OM4_025/WOA05_ptemp_salt_annual.v20141007.nc')
lat = WOA.latitude.rename({'latitude':'lat'})
lon = WOA.longitude.rename({'longitude':'lon'})
WOA_grid = xr.merge([lon,lat])


regridder = xe.Regridder(cesm_grid, WOA_grid, 'bilinear')

