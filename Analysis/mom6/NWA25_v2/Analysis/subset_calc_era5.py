# slice down the data
import xarray as xr
import os
import cftime
import numpy as np
from glob import glob
import os

import matplotlib as mpl
import matplotlib.pyplot as plt

# Functions for humidity borrowed and adapted from MetPy.calc: https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.html
def mixing_ratio(partial_press, total_press, molecular_weight_ratio=0.622):
    return (molecular_weight_ratio * partial_press
                / (total_press - partial_press))


def specific_humidity_from_mixing_ratio(mr):
    return mr / (1 + mr)


def saturation_vapor_pressure(temperature):
    sat_pressure_0c = 6.112e2 # Pa
    return sat_pressure_0c * np.exp(17.67 * (temperature - 273.15) # K -> C
                                        / (temperature - 29.65))   # K -> C

def saturation_mixing_ratio(total_press, temperature):
    return mixing_ratio(saturation_vapor_pressure(temperature), total_press)


era5_dict = {
            'ERA5_10m_u_component_of_wind':'u10',
            'ERA5_10m_v_component_of_wind':'v10',
            'ERA5_2m_temperature':'t2m',
            'ERA5_surface_solar_radiation_downwards':'ssrd',
            'ERA5_surface_thermal_radiation_downwards':'strd',
            'ERA5_total_rain_rate':'trr',
            'ERA5_mean_sea_level_pressure':'msl',
            'ERA5_2m_specific_humidity':'huss'
            }

df = xr.open_dataset('/okyanus/users/milicak/dataset/ERA5/global/ERA5_2m_temperature_2018.nc')
xinds = np.arange(1000,1300)
yinds = np.arange(120,350)
ds = df.isel(longitude=xinds,latitude=yinds)
dlon = np.copy(ds.longitude)
dlon[dlon>180] = dlon[dlon>180]-360

year = 1995
years = range(year,year+1)
#subset
era5dir = "/okyanus/users/milicak/dataset/ERA5/global/"
subdir = '/okyanus/users/milicak/dataset/ERA5/NWA25_v2/subset/'

for f in era5_dict.keys():
    print(f)
    for y in years:
        print(y)
        if f=='ERA5_total_rain_rate':
            #crr = xr.open_dataset(str(era5dir + 'ERA5_convective_rain_rate_' + str(y) + '.nc')).sel(latitude=slice(lat_min,lat_max), longitude=slice(lon_min,lon_max))
            #lsrr = xr.open_dataset(str(era5dir + 'ERA5_large_scale_rain_rate_' + str(y) + '.nc')).sel(latitude=slice(lat_min,lat_max), longitude=slice(lon_min,lon_max))
            crr = xr.open_dataset(str(era5dir + 'ERA5_convective_rain_rate_' + str(y) + '.nc')).isel(longitude=xinds,latitude=yinds)
            lsrr = xr.open_dataset(str(era5dir + 'ERA5_large_scale_rain_rate_' + str(y) + '.nc')).isel(longitude=xinds,latitude=yinds)
            #trr = crr.drop('crr')
            trr = xr.Dataset()
            trr['trr'] = crr['crr'] + lsrr['lsrr']
            trr['trr'].attrs = {'units': 'kg m-2 s-1','long_name': 'Total Rainfall Rate'}
            trr['longitude'] = dlon
            trr.trr.encoding = {k: v for k, v in crr.crr.encoding.items() if k in {'_FillValue', 'missing_value', 'dtype'}}
            all_vars = list(trr.data_vars.keys()) + list(trr.coords.keys())
            encodings = {v: {'_FillValue': 1.0e20} for v in all_vars}
            encodings['time'].update({'dtype':'float64', 'calendar': 'gregorian'})
            encodings['trr'].update({'dtype':'float32'})
            trr.to_netcdf(str(subdir + f + '_' + str(y) + ".nc"), format="NETCDF4_CLASSIC", encoding=encodings, unlimited_dims='time')
            crr.close()
            lsrr.close()
            trr.close()
        if f=='ERA5_2m_specific_humidity':
            #pair = xr.open_dataset(str(era5dir + 'ERA5_surface_pressure_' + str(y) + '.nc'))['sp'].sel(latitude=slice(lat_min,lat_max), longitude=slice(lon_min,lon_max)) # Pa
            #tdew = xr.open_dataset(str(era5dir + 'ERA5_2m_dewpoint_temperature_' + str(y) + '.nc'))['d2m'].sel(latitude=slice(lat_min,lat_max), longitude=slice(lon_min,lon_max)) # K
            pair = xr.open_dataset(str(era5dir + 'ERA5_surface_pressure_' + str(y) + '.nc'))['sp'].isel(longitude=xinds,latitude=yinds)
            tdew = xr.open_dataset(str(era5dir + 'ERA5_2m_dewpoint_temperature_' + str(y) + '.nc'))['d2m'].isel(longitude=xinds,latitude=yinds)

            smr = saturation_mixing_ratio(pair, tdew)
            sphum = specific_humidity_from_mixing_ratio(smr)

            sphum.name = 'huss'
            sphum = sphum.to_dataset()
            sphum['longitude'] = dlon

            # Remove all _FillValue
            all_vars = list(sphum.data_vars.keys()) + list(sphum.coords.keys())
            encodings = {v: {'_FillValue': None} for v in all_vars}

            # Also fix the time encoding
            encodings['time'].update({'dtype':'float64', 'calendar': 'gregorian', 'units': 'hours since 1900-01-01 00:00:00'})

            fout=str(subdir + f + '_' + str(y) + ".nc")
            sphum.to_netcdf(
                fout,
                format='NETCDF4_CLASSIC',
                engine='netcdf4',
                encoding=encodings,
                unlimited_dims=['time']
            )
            sphum.close()

        if 'total_rain_rate' not in f and 'specific_humidity' not in f:
            #ds=xr.open_dataset(str(era5dir + f + '_' + str(y) + ".nc")).sel(latitude=slice(lat_min,lat_max), longitude=slice(lon_min,lon_max))
            ds=xr.open_dataset(str(era5dir + f + '_' + str(y) + ".nc")).isel(longitude=xinds,latitude=yinds)
            ds['longitude'] = dlon
            ds.to_netcdf(str(subdir + f + '_' + str(y) + ".nc"),format="NETCDF4_CLASSIC")
            ds.close()
