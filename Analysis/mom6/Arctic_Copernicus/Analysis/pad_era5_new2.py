import xarray as xr
import os
import cftime
import numpy as np

datadir = '/okyanus/users/milicak/dataset/ERA5/Arctic_Copernicus/padded/'
outdir = '/okyanus/users/milicak/dataset/ERA5/Arctic_Copernicus/flooded/'

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

yr = 2002
for year in range(yr,yr+1):
    print(year)
    for f, f1 in era5_dict.items():
        print(f)
        # open the file for current year
        df = xr.open_dataset(f"{datadir}/{f}_{year}.nc")
        df = df.rename({'longitude':'x'})
        lon = np.zeros(df.x.shape[0]+1)
        lon[1:] = np.copy(df.x)
        lon[0] = -180.0
        ds = xr.Dataset({
            'longitude': xr.DataArray(
                        data   = lon,   # enter data here
                        dims   = ['longitude'],
                        coords = {'longitude': lon},
                        )})
        temp = np.zeros((df[f1].shape[0],df[f1].shape[1],df[f1].shape[2]+1))
        temp[:,:,1:] = np.copy(df[f1])
        temp[:,:,0] = np.copy(df[f1][:,:,-1])
        out = xr.Dataset({
            f1: xr.DataArray(
                        data   = temp,   # enter data here
                        dims   = ['time','latitude','longitude'],
                        coords = {'time': df.time, 'latitude': df.latitude,'longitude': ds.longitude},
                        )})
        all_vars = list(df.data_vars.keys()) + list(out.coords.keys())
        if f1=='ssrd' or f1=='strd':
            out[f1].attrs['units'] = 'W m-2'
        if f1=='huss':
            out[f1].attrs['dtype'] = 'float64'
            out[f1].attrs['standard_name'] = 'specific_humidity'
            out[f1].attrs['long_name'] = 'Near-Surface Specific Humidity'
            out[f1].attrs['coordinates'] = 'height'
            out[f1].attrs['units'] = '1'
            out['height'] = 2.0
            out['height'].attrs['units'] = "m"
            out['height'].attrs['axis'] = "Z"
            out['height'].attrs['positive'] = "up"
            out['height'].attrs['long_name'] = "height"
            out['height'].attrs['standard_name'] = "height"
        if f1=='t2m':
            out['height'] = 2.0
            out['height'].attrs['units'] = "m"
            out['height'].attrs['axis'] = "Z"
            out['height'].attrs['positive'] = "up"
            out['height'].attrs['long_name'] = "height"
            out['height'].attrs['standard_name'] = "height"
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
        out.to_netcdf(f'{outdir}{f}_{year}.nc', format="NETCDF4_CLASSIC", encoding=encodings, unlimited_dims='time')
        out.close()






