# slice down the data
import xarray as xr
import os
import cftime
import numpy as np
from glob import glob
import os

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

era5_dict = {'ERA5_10m_u_component_of_wind':'u10',
            'ERA5_10m_v_component_of_wind':'v10',
            'ERA5_2m_temperature':'t2m',
            'ERA5_surface_solar_radiation_downwards':'ssrd',
            'ERA5_surface_thermal_radiation_downwards':'strd',
            'ERA5_total_rain_rate':'trr',
            'ERA5_mean_sea_level_pressure':'msl',
            'ERA5_2m_specific_humidity':'huss'}

years=range(1996,1998)
#subset
era5dir = "/okyanus/users/milicak/dataset/ERA5/global/"
subdir = '/okyanus/users/milicak/dataset/ERA5/Arctic_Copernicus/subset/'
lon_min = 0;
lon_max = 360;
lat_min = 90;
lat_max = 25;

for f in era5_dict.keys():
    print(f)
    for y in years:
        print(y)
        if f=='ERA5_total_rain_rate':
            crr = xr.open_dataset(str(era5dir + 'ERA5_convective_rain_rate_' + str(y) + '.nc')).sel(latitude=slice(lat_min,lat_max),longitude=slice(lon_min,lon_max))
            lsrr = xr.open_dataset(str(era5dir + 'ERA5_large_scale_rain_rate_' + str(y) + '.nc')).sel(latitude=slice(lat_min,lat_max),longitude=slice(lon_min,lon_max))
            trr = crr.drop('crr')
            trr['trr'] = crr['crr'] + lsrr['lsrr']
            trr['trr'].attrs = {'units': 'kg m-2 s-1','long_name': 'Total Rainfall Rate'}
            trr.trr.encoding = {k: v for k, v in crr.crr.encoding.items() if k in {'_FillValue', 'missing_value', 'dtype'}}
            #trr.trr.encoding.update({'add_offset': None, 'scale_factor': None})
            trr.to_netcdf(str(subdir + f + '_' + str(y) + ".nc"), mode='w', format='NETCDF4_CLASSIC')
            crr.close()
            lsrr.close()
            trr.close()
        if f=='ERA5_2m_specific_humidity':
            pair = xr.open_dataset(str(era5dir + 'ERA5_surface_pressure_' + str(y) + '.nc'))['sp'].sel(latitude=slice(lat_min,lat_max),longitude=slice(lon_min,lon_max)) # Pa
            tdew = xr.open_dataset(str(era5dir + 'ERA5_2m_dewpoint_temperature_' + str(y) + '.nc'))['d2m'].sel(latitude=slice(lat_min,lat_max),longitude=slice(lon_min,lon_max)) # K

            smr = saturation_mixing_ratio(pair, tdew)
            sphum = specific_humidity_from_mixing_ratio(smr)

            sphum.name = 'huss'
            sphum = sphum.to_dataset()

            times = pair['time']
            latitudes = pair['latitude']
            longitudes = pair['longitude']
            # extend the longitude coordinate by 1
            datasets = []
            datasets.append(longitudes)
            # add 0.25 degrees to the final longitude value because the resolution of the ERA5 is 0.25 degrees
            datasets.append(longitudes[-1] + 0.25)
            longitudes = xr.concat(datasets, dim='longitude')
            f1 = 'huss'
            ds_ext1 = xr.Dataset({
                f1 : xr.DataArray(
                            data   = np.zeros((len(range(0,pair.dims['time'])),
                                               len(range(0,pair.dims['latitude'])),
                                               len(range(0,pair.dims['longitude']+1)))),   # enter data here
                            dims   = ['time', 'latitude', 'longitude'],
                            coords = {'time': times, 'latitude' : latitudes, 'longitude' : longitudes},
                            attrs  = {
                                'units'     : sphum[f1].attrs['units']
                                }
                            )})

            ds_ext1[f1].values[:,:,0:len(pair['longitude'])] = sphum[f1].values
            ds_ext1[f1].values[:,:,len(pair['longitude']):len(pair['longitude']) + 1] = sphum[f1][:,:,0:1].values

            # Remove all _FillValue
            all_vars = list(sphum.data_vars.keys()) + list(sphum.coords.keys())
            encodings = {v: {'_FillValue': None} for v in all_vars}

            # Also fix the time encoding
            encodings['time'].update({'dtype':'float64', 'calendar': 'gregorian', 'units': 'hours since 1900-01-01 00:00:00'})

            fout=str(subdir + f + '_' + str(y) + ".nc")
            # sphum.to_netcdf(
            ds_ext1.to_netcdf(
                fout,
                format='NETCDF4_CLASSIC',
                engine='netcdf4',
                encoding=encodings,
                unlimited_dims=['time']
            )
            sphum.close()
            ds_ext1.close()

        if 'total_rain_rate' not in f and 'specific_humidity' not in f:
            ds = xr.open_dataset(str(era5dir + f + '_' + str(y) + ".nc")).sel(latitude=slice(lat_min,lat_max),longitude=slice(lon_min,lon_max))
            times = ds['time']
            latitudes = ds['latitude']
            longitudes = ds['longitude']
            # extend the longitude coordinate by 1
            datasets = []
            datasets.append(longitudes)
            # add 0.25 degrees to the final longitude value because the resolution of the ERA5 is 0.25 degrees
            datasets.append(longitudes[-1] + 0.25)
            longitudes = xr.concat(datasets, dim='longitude')

            ds_ext = xr.Dataset({
                f : xr.DataArray(
                            data   = np.zeros((len(range(0,ds.dims['time'])),
                                               len(range(0,ds.dims['latitude'])),
                                               len(range(0,ds.dims['longitude']+1)))),   # enter data here
                            dims   = ['time', 'latitude', 'longitude'],
                            coords = {'time': times, 'latitude' : latitudes, 'longitude' : longitudes},
                            attrs  = {
                                'units'     : ds[f].attrs['units']
                                }
                            )})

            ds_ext[f].values[:,:,0:len(ds['longitude'])] = ds[f].values
            ds_ext[f].values[:,:,len(ds['longitude']):len(ds['longitude']) + 1] = ds[f][:,:,0:1].values

            ds.to_netcdf(str(subdir + f + '_' + str(y) + ".nc"),format="NETCDF4_CLASSIC")
            ds.close()

