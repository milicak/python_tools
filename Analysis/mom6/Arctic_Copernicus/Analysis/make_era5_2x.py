import numpy as np


root_dir = '/okyanus/users/milicak/dataset/ERA5/Arctic_Copernicus/flooded/'

year = 1996

ls1 = sorted(glob.glob(root_dir+'*'+str(year)+'*.nc'))

for fname in ls1:
    df = xr.open_dataset(fname)
    df = df.sel(time=slice("1995-12-31", "1996-02-01"))
    fout = fname[:-3]+'_2X.nc'
    print(fout)
    df.to_netcdf(fout,unlimited_dims=('time'))
