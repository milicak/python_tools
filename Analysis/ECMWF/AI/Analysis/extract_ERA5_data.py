import numpy as np

ls1 = sorted(glob.glob('/okyanus/users/milicak/dataset/ERA5/era5_global_2m_temperature_*.nc'))

for file in ls1:
    print(file)
    df = xr.open_dataset(file)
    min_lon = -1
    min_lat = 25.05
    max_lon = 60.05
    max_lat = 89.95

    mask_lon = (df.longitude >= min_lon) & (df.longitude <= max_lon)
    mask_lat = (df.latitude >= min_lat) & (df.latitude <= max_lat)

    cropped_ds1 = df.where(mask_lon & mask_lat, drop=True)

    min_lon = 280
    min_lat = 25.05
    max_lon = 360.5
    max_lat = 89.95

    mask_lon = (df.longitude >= min_lon) & (df.longitude <= max_lon)
    mask_lat = (df.latitude >= min_lat) & (df.latitude <= max_lat)

    cropped_ds2 = df.where(mask_lon & mask_lat, drop=True)
    cropped_ds2['longitude']=cropped_ds2['longitude']-360

    ds = xr.merge([cropped_ds2,cropped_ds1])
    dsdaily = ds.resample(time='1D').mean('time')
    fout = file[:67] + '_Europe.nc'
    dsdaily.to_netcdf(fout)
