import numpy as np

df = xr.open_dataset('/data/opa/bs-mod/upstream_bs-nrt/baseline/geodta/domain_cfg_MIv3.nc')

ds = df.isel(y=slice(17,261))
ds['jpjglo'] = np.int32(244)

ds['bathy_metry'][:,0:5,40:64] = 0
ds['top_level'][:,0:5,40:64] = 0
ds['bottom_level'][:,0:5,40:64] = 0

ds.attrs['DOMAIN_position_last'] = np.array([591,244],dtype=np.int32)
ds.attrs['DOMAIN_size_global'] = np.array([591,244],dtype=np.int32)
ds.attrs['DOMAIN_size_local'] = np.array([591,244],dtype=np.int32)

plt.pcolormesh(ds.bathy_metry[0,:,:].where(ds.bathy_metry[0,:,:]!=0));plt.colorbar();

ds.to_netcdf('/work/opa/mi19918/Projects/shared/domain_cfg_MIv4.nc')
