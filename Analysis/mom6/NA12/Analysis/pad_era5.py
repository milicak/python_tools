import xarray as xr
import os
import cftime
import numpy as np

datadir = '/okyanus/users/milicak/dataset/ERA5/NA12/subset/'
outdir = '/okyanus/users/milicak/dataset/ERA5/NA12/padded/'

era5_dict = {
            # 'ERA5_10m_u_component_of_wind':'u10',
            # 'ERA5_10m_v_component_of_wind':'v10',
            # 'ERA5_2m_temperature':'t2m',
            # 'ERA5_surface_solar_radiation_downwards':'ssrd',
            # 'ERA5_surface_thermal_radiation_downwards':'strd',
            'ERA5_total_rain_rate':'trr',
            # 'ERA5_mean_sea_level_pressure':'msl',
            # 'ERA5_2m_specific_humidity':'huss'
            }

for year in range(1995,2018):
    print(year)
    for f, f1 in era5_dict.items():
        print(f)
        # open the file for current year
        current = xr.open_dataset(f"{datadir}/{f}_{year}.nc")
        # next_data = xr.open_dataset(f"{datadir}/{f}_{year+1}.nc").isel(time=0)
        next_data = xr.open_dataset(f"{datadir}/{f}_{year+1}.nc").isel(time=slice(0,2))
        if year != 1995:
            previous = xr.open_dataset(f"{datadir}/{f}_{year-1}.nc").isel(time=-1)
            out = xr.concat([previous, current, next_data], dim="time")
            previous.close()
        else:
            out = xr.concat([current, next_data], dim="time")

        current.close()
        next_data.close()
        all_vars = list(out.data_vars.keys()) + list(out.coords.keys())
        encodings = {v: {'_FillValue': 1.0e20} for v in all_vars}
        encodings['time'].update({'dtype':'float64', 'calendar': 'gregorian'})
        out['time'].attrs['long_name'] = 'time'
        out['time'].attrs['standard_name'] = 'time'
        out['time'].attrs['axis'] = 'T'
        out['latitude'].attrs['long_name'] = 'Latitude'
        out['latitude'].attrs['units'] = 'degrees_north'
        out['latitude'].attrs['axis'] = 'Y'
        out['longitude'].attrs['long_name'] = 'Longitude'
        out['longitude'].attrs['units'] = 'degrees_east'
        out['longitude'].attrs['axis'] = 'X'
        out=out.transpose("time", "latitude", "longitude")
        # latitude needs to be reindexed for some reason
        out=out.reindex(latitude=list(reversed(out.latitude)))
        out.to_netcdf(f'{outdir}{f}_{year}.nc', format="NETCDF4_CLASSIC", encoding=encodings, unlimited_dims='time')
        out.close()






