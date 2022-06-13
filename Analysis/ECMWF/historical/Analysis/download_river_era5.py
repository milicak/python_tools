import cdsapi

c = cdsapi.Client()

c.retrieve(
    'cems-glofas-historical',
    {
        'variable': 'river_discharge_in_the_last_24_hours',
        'format': 'grib',
        'system_version': 'version_3_1',
        'hyear': '1995',
        'hmonth': [
            'april', 'august', 'december',
            'february', 'january', 'july',
            'june', 'march', 'may',
            'november', 'october', 'september',
        ],
        'hydrological_model': 'lisflood',
        'hday': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
            '13', '14', '15',
            '16', '17', '18',
            '19', '20', '21',
            '22', '23', '24',
            '25', '26', '27',
            '28', '29', '30',
            '31',
        ],
        # 'area': [89.95, -100, -35,40,],
        'product_type': 'consolidated',
    },
    'glofas-era5_1995.grib')

