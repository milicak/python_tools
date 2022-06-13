import os

variables = [
             '10m_u_component_of_wind','10m_v_component_of_wind','2m_dewpoint_temperature',
             '2m_temperature','mean_sea_level_pressure','snowfall',
             'surface_pressure','surface_solar_radiation_downwards','surface_thermal_radiation_downwards',
             'convective_rain_rate','large_scale_rain_rate',
             'total_cloud_cover','total_precipitation'
         ]

variables = [
             'snowfall'
             'surface_pressure',
             'surface_solar_radiation_downwards',
             'surface_thermal_radiation_downwards',
             'convective_rain_rate','large_scale_rain_rate',
             'total_cloud_cover','total_precipitation'
         ]
         #
root_folder  = '/okyanus/users/milicak/dataset/ERA5/global/'

# for year in range(1997,2003):
years = 2018
for year in range(years,years+1):
    for var in variables:
        print(var)
        cmnd = 'scp milicak@antares.esm.rutgers.edu:/Volumes/P5/ERA5/*' + var + '_' + str(year) + '.nc ' + root_folder + '.'
        os.system(cmnd)
