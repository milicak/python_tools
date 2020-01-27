import cdsapi

c = cdsapi.Client()
c.retrieve('reanalysis-era5-pressure-levels', {
        'variable'      : '2m temperature',
        'levtype'       : 'sfc',
        'product_type'  : 'reanalysis',
        'year'          : '2008',
        'month'         : '01',
        'day'           : '01',
        'time'          : '12:00',
        'format'        : 'netcdf' # Supported format: grib and netcdf. Default: grib
    }, 'test.nc')
