import cdsapi

c = cdsapi.Client()

root_folder = '/okyanus/users/milicak/dataset/ERA5/'

variables = ['10m_u_component_of_wind','10m_v_component_of_wind','2m_dewpoint_temperature',
'2m_temperature','evaporation','mean_sea_level_pressure','surface_latent_heat_flux',
'surface_net_solar_radiation','surface_net_thermal_radiation','surface_sensible_heat_flux',
'surface_thermal_radiation_downwards','total_cloud_cover','total_precipitation']

# variables = ['2m_dewpoint_temperature']
# variables = ['evaporation']
# variables = ['mean_sea_level_pressure']
# variables = ['surface_latent_heat_flux']
# variables = ['surface_net_solar_radiation']
# variables = ['surface_net_thermal_radiation']
# variables = ['surface_sensible_heat_flux']
# variables = ['total_cloud_cover']
# variables = ['snowfall']

# For MITgcm exf package
# variables = ['10m_u_component_of_wind']
# variables = ['10m_v_component_of_wind']
# variables = ['2m_temperature']
# variables = ['surface_thermal_radiation_downwards']
# variables = ['surface_solar_radiation_downwards']
# variables = ['total_precipitation']
# variables = ['surface_runoff']
variables = ['mean_surface_runoff_rate']

for year in range(2019, 2020):
# for year in range(2003, 2004):
# for year in range(2004, 2019):
    for var in variables:
        filename = root_folder + 'era5_global_' + var + '_' + str(year) + '.nc'
        c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type':'reanalysis',
                'format':'netcdf',
                'variable':[
                    var
                ],
                # 'area'    : "48.50/25.00/37.00/43.00",
                'year':[
                    year
                ],
                'month':[
                    '01','02','03',
                    '04','05','06',
                    '07','08','09',
                    '10','11','12'
                ],
                'day':[
                    '01','02','03',
                    '04','05','06',
                    '07','08','09',
                    '10','11','12',
                    '13','14','15',
                    '16','17','18',
                    '19','20','21',
                    '22','23','24',
                    '25','26','27',
                    '28','29','30',
                    '31'
                ],
                # 'time':[
                #     '00:00','01:00','02:00',
                #     '03:00','04:00','05:00',
                #     '06:00','07:00','08:00',
                #     '09:00','10:00','11:00',
                #     '12:00','13:00','14:00',
                #     '15:00','16:00','17:00',
                #     '18:00','19:00','20:00',
                #     '21:00','22:00','23:00',
                # ]
                'time':[
                    '00:00',
                    '06:00',
                    '12:00',
                    '18:00',
                ]
            },
            filename)




# '10m_u_component_of_wind','10m_v_component_of_wind','2m_dewpoint_temperature',
# '2m_temperature','evaporation','mean_sea_level_pressure','surface_latent_heat_flux',
# 'surface_net_solar_radiation','surface_net_thermal_radiation','surface_sensible_heat_flux',
# 'surface_thermal_radiation_downwards','total_cloud_cover','total_precipitation'
