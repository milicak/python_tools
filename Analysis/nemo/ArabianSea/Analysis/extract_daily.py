import numpy as np
import glob
import xarray as xr
import os
import cftime
import numpy as np

datadir = '/okyanus/users/milicak/dataset/ERA5/ArabianSea/padded/'
outdir = '/okyanus/users/milicak/dataset/ERA5/ArabianSea/daily/'


dailyval = np.int16(24)
yr = 1993
# for year in range(yr,yr+1):

for year in range(yr,2021):
    print(year)
    ls1=sorted(glob.glob(datadir+'*'+str(year)+'*'))
    df = xr.open_mfdataset(ls1)
    for days in np.arange(1,df.time.dt.dayofyear[-1]+1):
        out = df.isel(time=(df.time.dt.dayofyear == days))
        fout = outdir + 'era5_y' + str(year) + 'm' + str(np.copy(out.time.dt.month[0])).zfill(2) + 'd' + str(np.copy(out.time.dt.day[0])).zfill(2) + '.nc'
        print(fout)
        out = out.rename({'time': 'time_counter'})
        all_vars = list(out.data_vars.keys()) + list(out.coords.keys())
        encodings = {v: {'_FillValue': None} for v in all_vars}
        encodings['time_counter'].update({'dtype':'float64', 'calendar': 'gregorian'})
        out.to_netcdf(
                fout,
                format='NETCDF4_CLASSIC',
                engine='netcdf4',
                encoding=encodings,
                unlimited_dims=['time_counter'])
        out.close()
        out.to_netcdf(fout, format="NETCDF4_CLASSIC", encoding=encodings,
                      unlimited_dims='time_counter')
        out.close()



