import xesmf as xe
import ESMF

df = xr.open_dataset('/work/opa/mi19918/Projects/uTSS_SHYFEM/work/monthly/uTSS_lobc_chunk_monthly_2020_2021.nos.nc')
ds = xr.open_dataset('NEMO_BS_boundary_coords.nc')
# dss=xr.open_dataset('/data/opa/bs-mod/upstream_bs-nrt/bdydta_refresh_utss2021/bdyT_tra_m07.nc')

#section 1
ds_locs = xr.Dataset()
ds_locs["lon"] = xr.DataArray(
    data=np.copy(ds.nav_lon[0,0:8]), dims=("locations1")
)
ds_locs["lat"] = xr.DataArray(data=np.copy(ds.nav_lat[0,0:8]),
                              dims=("locations2"))

dd = df.isel(month=0)
dd = dd.isel(level=0)
regridder = xe.Regridder(dd, ds_locs, 'nearest_s2d',locstream_in=True)

temp = np.zeros((12,92,8))
salt = np.zeros((12,92,8))
for time in range(0,12):
    print(time)
    for kind in range(0,92):
        print(kind)
        dd = df.isel(month=time)
        dd = dd.isel(level=kind)
        # regridder = xe.Regridder(dd, ds_locs, 'nearest_s2d',locstream_in=True)
        df_locs = regridder(dd)
        temp[time,kind,:] = np.copy(df_locs.temperature[:,0])
        salt[time,kind,:] = np.copy(df_locs.salinity[:,0])

temp[temp==0]=np.nan
salt[salt==0]=np.nan

# put data into a dataset
dsf = xr.Dataset({
    'temperature': xr.DataArray(
                data   = temp,
                dims   = ['time','depth','locations'],
        coords = {'time': np.copy(df.month), 'depth': np.copy(df.level), 'locations':
                  np.copy(ds_locs.locations1)},
                attrs  = {
                    # '_FillValue': -999.9,
                    'units'     : 'Celcius'
                    }
                ),
    'salinity': xr.DataArray(
                data   = salt,
                dims   = ['time','depth','locations'],
        coords = {'time': np.copy(df.month), 'depth': np.copy(df.level), 'locations':
                  np.copy(ds_locs.locations1)},
                attrs  = {
                    # '_FillValue': -999.9,
                    'units'     : 'psu'
                    }
                )
            },
        # attrs = {'example_attr': 'this is a global attribute'}
    )

