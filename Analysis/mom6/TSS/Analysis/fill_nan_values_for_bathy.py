import xarray as xr
import scipy.io

# firs run make_inicond_shyfem_TSssh_interpolated.py

mom_dir = '/okyanus/users/milicak/dataset/MOM6/TSS/'
path_regional_grid = mom_dir + './ocean_hgrid.nc'
df2 = xr.open_dataset(path_regional_grid)
lon_mom = np.array(np.copy(df2['x'][1::2,1::2]), dtype=float)
lat_mom = np.array(np.copy(df2['y'][1::2,1::2]), dtype=float)

df = xr.open_dataset('mom6_TSS_salinity.nc')
df['yh'] = lon_mom[:,0]
# df['xh'] = lat_mom[0,:]
df['xh'] = np.flipud(lat_mom[0,:])
dnm = df.salinity.interpolate_na(dim="yh", method="nearest",fill_value="extrapolate")
dnm2 = dnm.interpolate_na(dim="xh", method="nearest",fill_value="extrapolate")
dnm2 = dnm2.expand_dims(dim={"time": [np.datetime64('1996-01-01')]}, axis=0)
dnm2 = dnm2[:,:,:,:-1]
dnm2 = dnm2.transpose("time","zl","yh","xh")

df = xr.open_dataset('mom6_TSS_temperature.nc')
df['yh'] = lon_mom[:,0]
df['xh'] = np.flipud(lat_mom[0,:])
dnm = df.temperature.interpolate_na(dim="yh", method="nearest",fill_value="extrapolate")
dnm3 = dnm.interpolate_na(dim="xh", method="nearest",fill_value="extrapolate")
dnm3 = dnm3.expand_dims(dim={"time": [np.datetime64('1996-01-01')]}, axis=0)
dnm3 = dnm3[:,:,:,:-1]
dnm3 = dnm3.transpose("time","zl","yh","xh")


df1 = dnm2.to_dataset(name='salt')
df1['temp'] = dnm3
df1['salt'] = df1.salt+0.15
all_vars = list(df1.data_vars.keys()) + list(df1.coords.keys())
encodings = {v: {'_FillValue': 1.0e20} for v in all_vars}
encodings['time'].update({'dtype':'float64', 'calendar': 'gregorian'})
df1['time'].attrs['long_name'] = 'time'
df1['time'].attrs['standard_name'] = 'time'
df1['time'].attrs['axis'] = 'T'

df1['zl'].attrs['long_name'] = 'Layer pseudo-depth, -z*'
df1['zl'].attrs['units'] = 'meter'
df1['zl'].attrs['cartesian_axis'] = 'Z'

df1.to_netcdf('tempsalt_IC.nc',format="NETCDF4_CLASSIC", encoding=encodings, unlimited_dims='time')

# Create a mosaic file
t, nk, nj, ni = dnm2.shape
time = 17.5
fout  = mom_dir + 'tempsalt_IC.nc'
rg = scipy.io.netcdf_file(fout,'w')
# Dimensions
# rg.createDimension('time', None)
rg.createDimension('time', 1)
rg.createDimension('depth',nk)
rg.createDimension('longitude',ni)
rg.createDimension('latitude',nj)
# Variables
hnx = rg.createVariable('longitude', 'float32', ('longitude',))
hnx.units = 'degrees east'
hnx.standard_name = 'longitude'
hny = rg.createVariable('latitude', 'float32', ('latitude',))
hny.units = 'degrees north'
hny.standard_name = 'latitude'
hz = rg.createVariable('depth','float32',('depth',))
hz.units = 'meters'
hz._CoordinateZisPositive = 'down'
tempvar  = rg.createVariable('temp','float32',('time','depth','latitude','longitude',))
tempvar.units = 'celcius'
tempvar.missing_val = 1e20
tempvar._FillValue = 1e20
saltvar  = rg.createVariable('salt','float32',('time','depth','latitude','longitude',))
saltvar.units = 'psu'
saltvar.missing_val = 1e20
saltvar._FillValue = 1e20
htime = rg.createVariable('time','float32',('time',))
# htime = rg.createVariable('time', 'int', ('time'))
htime.units = 'days since 1996-01-01 00:00:00'
# Values
# hnx[:] = np.copy(lon_mom[0,:])
# hny[:] = np.copy(lat_mom[:,0])
# hz[:] = np.copy(df.zl[:-1])
tempvar[:] = np.copy(dnm3)
saltvar[:] = np.copy(dnm2)
htime = time
rg.close()

