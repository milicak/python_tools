import xarray as xr
import os
import cftime
import numpy as np
import glob

datadir = '/okyanus/users/milicak/dataset/ERA5/BlackSea_nemo/subset/'
outdir = '/okyanus/users/milicak/dataset/ERA5/BlackSea_nemo/padded/'

era5_dict = {
            'ERA5_10m_u_component_of_wind':'u10',
            'ERA5_10m_v_component_of_wind':'v10',
            'ERA5_2m_temperature':'t2m',
            'ERA5_surface_solar_radiation_downwards':'ssrd',
            'ERA5_surface_thermal_radiation_downwards':'strd',
            'ERA5_total_rain_rate':'trr',
            'ERA5_mean_sea_level_pressure':'msl',
            'ERA5_2m_specific_humidity':'huss',
            'ERA5_snowfall':'sf',
            }

for year in range(1992,2020):
    print(year)
    # open the file for current year
    ls1 = sorted(glob.glob(datadir+'*'+str(year)+'*.nc'))
    out2 = xr.open_mfdataset(ls1)
    # out2 = xr.open_dataset(f"{datadir}/{f}_{year}.nc")
    for days in np.arange(1,out2.time.dt.dayofyear[-1]+1):
        out = out2.isel(time=(out2.time.dt.dayofyear == days))
        for f, f1 in era5_dict.items():
            if f1=='ssrd' or f1=='strd':
                # convert radiation from J/m2 to W/m2: https://confluence.ecmwf.int/pages/viewpage.action?pageId=155337784
                out[f1].values = out[f1].values/3600.0
                out[f1].attrs['units'] = 'W m-2'
            if f1=='huss':
                out[f1].attrs['dtype'] = 'float64'
                out[f1].attrs['standard_name'] = 'specific_humidity'
                out[f1].attrs['units'] = 'kg/kg'

        # latitude needs to be reindexed for some reason
        out = out.reindex(latitude=list(reversed(out.latitude)))
        fout = outdir + 'era5_y' + str(year) + 'm' + str(np.copy(out.time.dt.month[0])).zfill(2) + 'd' + str(np.copy(out.time.dt.day[0])).zfill(2) + '.nc'
        print(fout)
        out = out.rename({'time': 'time_counter'})
        all_vars = list(out.data_vars.keys()) + list(out.coords.keys())
        encodings = {v: {'_FillValue': None} for v in all_vars}
        encodings['time_counter'].update({'dtype':'float64', 'calendar': 'gregorian'})
        out.to_netcdf(
                fout,
                format='NETCDF4_CLASSIC',
                engine='netcdf4',
                encoding=encodings,
                unlimited_dims=['time_counter'])
        out.close()






