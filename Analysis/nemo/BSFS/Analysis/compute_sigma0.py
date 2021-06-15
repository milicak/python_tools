import numpy as np
import xarray as xr
import gsw

root_folder = '/okyanus/users/milicak/dataset/CMCC/BS/reanalysis/month/'

# for year in np.arange(1992,2019):
for year in np.arange(2018,2019):
    path_salt = root_folder + np.str(year) + '/*PSAL-BSe2r2*'
    path_temp = root_folder + np.str(year) + '/*TEMP-BSe2r2*'
    list1 = sorted(glob.glob(path_salt))
    list2 = sorted(glob.glob(path_temp))
    for month in np.arange(0,12):
        print(month, year)
        salt = xr.open_mfdataset(list1[month])
        temp = xr.open_mfdataset(list2[month])
        sigma = gsw.sigma0(salt.vosaline,temp.votemper)
        filename = list1[month]
        newfile = filename[:78]+'SIGM'+filename[82:]
        df = xr.DataArray(sigma, dims=["time", "deptht", "lat", "lon"],
                         coords=[salt.time, salt.depth, salt.lat, salt.lon])
        df = df.to_dataset(name='sigma')
        df.to_netcdf(newfile)




