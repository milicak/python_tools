import numpy as np
import glob
from datetime import timedelta

root_folder = '/okyanus/users/milicak/dataset/SODA3_12_21/sodafiles/'
year = 1980

ls1 = sorted(glob.glob(root_folder+'soda3.12.2_5dy_ocean_reg_'+str(year)+'*'))
ls2 = sorted(glob.glob(root_folder+'soda3.12.2_5dy_ocean_reg_'+str(year-1)+'*'))
ls3 = sorted(glob.glob(root_folder+'soda3.12.2_5dy_ocean_reg_'+str(year+1)+'*'))


if year == 1980:
    lsy = ls1+[ls3[0]]
else:
    lsy = [ls2[-1]]+ls1+[ls3[0]]

df = xr.open_mfdataset(lsy)


variables = [
    'temp',
    'salt',
    'u',
    'v',
    'ssh',
]
df = df[variables]

ds = df.resample(time='1m').mean(dim='time')
delta = timedelta(days=-15)
delta = timedelta(days=-30)
ds = df.resample(time='1m',loffset=delta).mean('time')

# all_vars = list(ds.data_vars.keys()) + list(ds.coords.keys())
# encodings = {v: {'_FillValue': None} for v in all_vars}
# encodings = {v: {'_FillValue': 1e20} for v in all_vars}

encodings = {'xt_ocean': {'zlib': False, '_FillValue': False},
            'yt_ocean': {'zlib': False, '_FillValue': False},
            'st_ocean': {'zlib': False, '_FillValue': False},
            'time': {'zlib': False, '_FillValue': False},
            'temp': {'_FillValue': 1e20},
            'salt': {'_FillValue': 1e20},
            'u': {'_FillValue': 1e20},
            'v': {'_FillValue': 1e20},
            'ssh': {'_FillValue': 1e20},
            }

# Make sure time has the right units and datatype
# otherwise it will become an int and MOM will fail.
encodings['time'].update({
        'units': 'days since 1950-01-01',
        'dtype': np.float,
        'calendar': 'NOLEAP'})

# Write out
out_file = '/okyanus/users/milicak/dataset/SODA3_12_21/SODA_monthly_sponge_' + str(year) + '.nc'

ds.to_netcdf(
    out_file,
    unlimited_dims=['time'],
    # format='NETCDF3_64BIT',
    encoding=encodings,
    engine='netcdf4')

ds.close()

# run following commands
# ncatted -a cartesian_axis,xt_ocean,c,c,"X" ~/dataset/SODA3_12_21/SODA_monthly_sponge_1980.nc
# ncatted -a cartesian_axis,yt_ocean,c,c,"Y" ~/dataset/SODA3_12_21/SODA_monthly_sponge_1980.nc
# ncatted -a cartesian_axis,st_ocean,c,c,"Z" ~/dataset/SODA3_12_21/SODA_monthly_sponge_1980.nc


