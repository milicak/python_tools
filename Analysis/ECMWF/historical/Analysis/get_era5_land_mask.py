import cdsapi

c = cdsapi.Client()

c.retrieve(
'reanalysis-era5-land',
{
'format': 'netcdf',
'variable': [
'land_sea_mask', 'orography'
],
'year': '1981',
'month': '01',
'day': '01',
'time': '01:00',
},
'download.nc')
