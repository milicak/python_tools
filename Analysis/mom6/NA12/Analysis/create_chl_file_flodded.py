import numpy as np
# import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import xarray as xr
from HCtFlood import kara as flood
import cftime

df = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/NA12/seawifs-clim-1997-2010.smoothed.nc')

flooded = flood.flood_kara(df['chlor_a'], xdim='i', ydim='j').drop('z').squeeze()
flooded.isel(time=0).plot()

hgrid = xr.open_dataset("/okyanus/users/milicak/dataset/MOM6/NA12/ocean_hgrid.nc")
# create a translator so xesmf knows the names of these variables in our hgrid
mom_grid = {
    'lon': hgrid.x[1::2, 1::2],
    'lon_b': hgrid.x[::2, ::2],
    'lat': hgrid.y[1::2, 1::2],
    'lat_b': hgrid.y[::2, ::2]
}
mom_area = (hgrid.area[::2, ::2] + hgrid.area[1::2, 1::2]) + (hgrid.area[1::2, ::2] + hgrid.area[::2, 1::2])

# reorganize our dimensions
ds = flooded.to_dataset(name='chlor_a')
ds = ds.rename({'i': 'xh', 'j': 'yh'})

# add in variables for our final file
ds['area'] = (('yh', 'xh'), mom_area.values)
ds['lon'] = (('yh', 'xh'), df['lon'].values)
ds['lat'] = (('yh', 'xh'), df['lat'].values)
ds['lon_crnr'] = (('yq', 'xq'), mom_grid['lon_b'].values)
ds['lat_crnr'] = (('yq', 'xq'), mom_grid['lat_b'].values)
ds['xh'] = (('xh'), np.arange(len(ds['xh']), dtype=np.int32))
ds['yh'] = (('yh'), np.arange(len(ds['yh']), dtype=np.int32))
ds['xh'].attrs['cartesian_axis'] = 'X'
ds['yh'].attrs['cartesian_axis'] = 'Y'


# Reassign time values. MOM6 is expecting a float 32 gregorian time dimension, therefore we use xarray cftime. Why the 732 hours? That's 30.5 days and roughly gives us the middle of each month
ds = ds.assign_coords(time=("time", xr.cftime_range(start=cftime.DatetimeNoLeap(1, 1, 15, 12, 0, 0, 0),end=cftime.DatetimeNoLeap(1, 12, 31, 23, 0, 0, 0), freq="732H")))
ds.time.attrs['modulo'] = ' '
ds.time.attrs['long_name'] = 'Time'
ds.time.attrs['cartesian_axis'] = 'T'
# ds.time.attrs['units'] = 'days since 0001-01-15'

# update the fill value for all variables attributes
all_vars = list(ds.data_vars.keys()) + list(ds.coords.keys())
encodings = {v: {'_FillValue': None} for v in all_vars}
encodings['time'].update({'dtype':'float64'})
encodings['time'].update({'units':'days since 0001-01-15'})

# ensure the data type is float 32 for the following variables
ds['lon'] = ds['lon'].astype(dtype='float32')
ds['lat'] = ds['lat'].astype(dtype='float32')
ds['area'] = ds['area'].astype(dtype='float32')
ds['lon_crnr'] = ds['lon_crnr'].astype(dtype='float32')
ds['lat_crnr'] = ds['lat_crnr'].astype(dtype='float32')

# save the file for MOM6. You can also save as `chl_mom6.nc`, as that is the default file name for this. Note that these parameters are controlled in MOM6 opacity parameters in MOM_input.
# This file MUST be output as netcdf3! MOM6 does not like netcdf4 (as of December 2021). This file should be placed in your INPUT folder (this is the default location where MOM6 searches for the file)

ds.to_netcdf('/okyanus/users/milicak/dataset/MOM6/NA12/seawifs-clim.1997.2010.nwa12_v2.nc',
    encoding=encodings,
    format='NETCDF3_64BIT',
    unlimited_dims='time'
            )

