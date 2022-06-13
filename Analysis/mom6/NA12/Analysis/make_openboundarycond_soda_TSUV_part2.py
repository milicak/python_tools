import numpy as np
import glob
import datetime

# ls1=sorted(glob.glob('/okyanus/users/milicak/dataset/SODA3_12_21/NA12_OBC/*1605*'))

root_folder = '/okyanus/users/milicak/dataset/SODA3_12_21/NA12_OBC/'

variables = ['ssh_north1605_obc','ssh_south1605_obc',
             'tracers_north_1605_obc','tracers_south_1605_obc',
             'uv_north1605_obc','uv_south1605_obc']


attrs = {'units': 'days since 1980-01-01','calendar': 'gregorian',
         'dtype':'float64'}

for var in variables:
    fname = root_folder + var + '.nc'
    df = xr.open_dataset(fname)
    # out3 = df.isel(time=-1)
    out3 = df
    # newtime = xr.Dataset({'time': datetime.datetime(2017, 12, 24)})
    ds = xr.Dataset({'time': ('time', [13872], attrs)})
    newtime=xr.decode_cf(ds)
    out3['time'] = newtime.time
    outname = root_folder + var[:-5]+'6_obc.nc'
    print(outname)
    # out3.time.encoding['calendar']='gregorian'
    # out3.time.encoding['units'] = 'days since 1980-01-01'
    out3.to_netcdf(outname,unlimited_dims='time',format='NETCDF3_CLASSIC')
    # newtime = xr.Dataset({'time': datetime.datetime(2017, 12, 29)})
    ds = xr.Dataset({'time': ('time', [13877], attrs)})
    newtime=xr.decode_cf(ds)
    out3['time'] = newtime.time
    # out3.time.encoding['calendar']='gregorian'
    # out3.time.encoding['units'] = 'days since 1980-01-01'
    outname = root_folder + var[:-5]+'7_obc.nc'
    print(outname)
    out3.to_netcdf(outname,unlimited_dims='time',format='NETCDF3_CLASSIC')
    # newtime = xr.Dataset({'time': datetime.datetime(2018, 1, 3)})
    # ds = xr.Dataset({'time': ('time', [13872, 13877, 13882], attrs)})
    ds = xr.Dataset({'time': ('time', [13882], attrs)})
    newtime=xr.decode_cf(ds)
    out3['time'] = newtime.time
    # out3.time.encoding['calendar']='gregorian'
    # out3.time.encoding['units'] = 'days since 1980-01-01'
    outname = root_folder + var[:-5]+'8_obc.nc'
    print(outname)
    out3.to_netcdf(outname,unlimited_dims='time',format='NETCDF3_CLASSIC')
