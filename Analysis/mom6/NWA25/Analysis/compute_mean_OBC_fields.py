import numpy as np
import xarray as xr
root_folder = '/Volumes/A1/workdir/milicak/datasets/MOM6/NWA25/OBC/'

ls1 = sorted(glob.glob(root_folder+'*.nc'))

for idx,fname in enumerate(ls1):
    ds = xr.open_dataset(fname)
    if "east" in fname:
        for varname in list(ds.keys()):
            tmp = np.copy(ds[varname].mean('time'))
            if(len(tmp.shape)==2):
                dnm = np.tile(tmp,(368,1,1))
                dstmp = xr.DataArray(dnm, coords={'time': ds.time,
                    'locations': ds.locations, 'dim_0': ds.dim_0,},
                     dims=['time', 'locations', 'dim_0'])
            else:
                dnm = np.tile(tmp,(368,1,1,1))
                dstmp = xr.DataArray(dnm, coords={'time': ds.time, 'z_l': ds.z_l, 
                    'locations': ds.locations, 'dim_0': ds.dim_0,},
                     dims=['time', 'z_l', 'locations', 'dim_0'])

            print('mehmet',varname)
            ds[varname] = dstmp
    else:
        for varname in list(ds.keys()):
            tmp = np.copy(ds[varname].mean('time'))
            if(len(tmp.shape)==2):
                dnm = np.tile(tmp,(368,1,1))
                dstmp = xr.DataArray(dnm, coords={'time': ds.time,
                    'dim_0': ds.dim_0, 'locations': ds.locations},
                     dims=['time', 'dim_0', 'locations'])
            else:
                dnm = np.tile(tmp,(368,1,1,1))
                dstmp = xr.DataArray(dnm, coords={'time': ds.time, 'z_l': ds.z_l, 
                    'dim_0': ds.dim_0, 'locations': ds.locations},
                     dims=['time', 'z_l', 'dim_0', 'locations'])

            print('ilicak',varname)
            ds[varname] = dstmp

    for v in ds:
        ds[v].encoding['_FillValue']=1.e20
    
    fout = fname[:-3]+'_mean.nc'
    ds['dim_0']=ds.dim_0.astype('int32')
    ds['locations']=ds.locations.astype('int32')
    print(fout)
    ds.to_netcdf(fout,unlimited_dims=('time'))
    ds.close()

