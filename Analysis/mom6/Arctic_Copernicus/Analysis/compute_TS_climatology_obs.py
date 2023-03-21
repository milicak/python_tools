import numpy as np
import xarray as xr
import xesmf as xe


root_folder = '/archive/milicak/MOM6-examples/Projects/Arctic_Copernicus/'

gr = xr.open_dataset(root_folder + 'ocean_geometry.nc')

momNA12= xr.Dataset()
momNA12['lon'] = gr['geolon']
momNA12['lat'] = gr['geolat']

dft = xr.open_dataset('/archive/milicak/dataset/WOA13/woa13_decav_t00_04v2.nc',decode_times=False)
dfs = xr.open_dataset('/archive/milicak/dataset/WOA13/woa13_decav_s00_04v2.nc',decode_times=False)
dfwoa= xr.Dataset()
dfwoa['lon'] = dfs['lon']
dfwoa['lat'] = dfs['lat']


# Calculate remapping weights
# Using nearest neighbor - other options could be used here , e.g. bilinear.
regrid_woa = xe.Regridder(dfwoa, momNA12, 'patch',
                       periodic=False,
                       filename='/archive/milicak/MOM6-examples/Projects/Arctic_Copernicus/regrid_woa_clim.nc',
                          reuse_weights=True)


temp_woa = regrid_woa(dft['t_an'][0,:,:,:])
dnm = temp_woa
dnm = dnm.fillna(1e3)
dnm = dnm.where(dnm!=0)
dnm = dnm.bfill(('lath'))
dnm = dnm.where(dnm!=1e3)
temp_woa = dnm

salt_woa = regrid_woa(dfs['s_an'][0,:,:,:])
dnm = salt_woa
dnm = dnm.fillna(1e3)
dnm = dnm.where(dnm!=0)
dnm = dnm.bfill(('lath'))
dnm = dnm.where(dnm!=1e3)
salt_woa = dnm

df_woa = temp_woa.to_dataset(name='temp_woa')
ds_woa = salt_woa.to_dataset(name='salt_woa')

df_woa['salt_woa'] = ds_woa['salt_woa']
df_woa.to_netcdf(root_folder + 'TS_woa_Arctic_Copernicus.nc')



