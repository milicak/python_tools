import cdsapi

c = cdsapi.Client()
c.retrieve('reanalysis-era5-complete', {
    'class': 'ea',
    'date': '2016-01-01',
    'expver': '1',
    'levtype': 'sfc',
    'param': '167.128',
    'step': '0',
    'stream': 'oper',
    'time': '06:00:00/18:00:00',
    'type': 'an',
    # 'format': 'netcdf', # Supported format: grib and netcdf. Default: grib
}, 'output.grb')
