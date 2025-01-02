
root_folder = '/okyanus/users/milicak//dataset/NEMO/TSS_BS_AGRIF/'
df1 = xr.open_dataset(root_folder + 'TS_climatology_BlackSea.nc')
df2 = xr.open_dataset(root_folder + 'TS_climatology_Marmara.nc')
df3 = xr.open_dataset(root_folder + 'TSy1993_m01_d01.nc')

# overwrite marmara sea kocaeliregion


path_regional_grid = root_folder + './domain_cfg.nc'
gr = xr.open_dataset(path_regional_grid)

df3['time_counter']  = df1['time_counter']


df1['vosaline'].loc[dict(x=slice(0,231))] = df3['so'].isel(x=slice(0,231))
df1['thetao'].loc[dict(x=slice(0,231))] = df3['thetao'].isel(x=slice(0,231))

df1['vosaline'].loc[dict(x=slice(231,365),y=slice(74,118))]=df2['vosaline'].isel(x=slice(231,365),y=slice(74,118))
df1['thetao'].loc[dict(x=slice(231,365),y=slice(74,118))]=df2['thetao'].isel(x=slice(231,365),y=slice(74,118))

tmp = df1

all_vars = list(tmp.data_vars.keys()) + list(tmp.coords.keys())
encodings = {v: {'_FillValue': None} for v in all_vars}
encodings['time_counter'].update({'dtype':'float64', 'calendar': 'gregorian'})
ftmp = root_folder + 'TSmergedy1993_m01_d01.nc'
tmp.to_netcdf(
        ftmp,
        format='NETCDF4_CLASSIC',
        engine='netcdf4',
        unlimited_dims=['time_counter'])
tmp.close()

