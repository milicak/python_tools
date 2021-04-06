import xarray as xr
import numpy as np

root_folder = '/okyanus/users/milicak/dataset/ERA5/'

varname = 'era5_global_total_precipitation_24h_'
var = 'tp'
var2 = 'rain'
coeff = 3600.0
# varname = 'era5_global_surface_thermal_radiation_downwards_'
# var = 'strd'
# var2 = 'dlw'
# coeff = 3600.0
# varname = 'era5_global_surface_solar_radiation_downwards_'
# var = 'ssrd'
# var2 = 'dsw'
# coeff = 3600.0

for year in range(2004,2019):
    fname = root_folder + varname + str(year) + '.nc'
    df = xr.open_dataset(fname)
    ds = df.resample(time='6H')
    dsm = ds.mean('time')
    dnm = np.flip(np.copy(dsm[var])/coeff,axis=1)
    outputfile = root_folder + 'ERA5_' + var2 + '_' + str(year)
    print(outputfile)
    dnm.astype('>f4').tofile(outputfile)




