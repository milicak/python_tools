from vcr import utils, conserve
from HCtFlood.kara import flood_kara
import xesmf as xe
# import GFDL_xr


def convert_lon(ds, lon='xh'):
    ds.coords[lon] = (ds.coords[lon] + 180) % 360 - 180
    ds = ds.sortby(ds[lon])
    return ds


ds = xr.open_dataset('/archive/milicak/MOM6-examples/Projects/patchy_NA/work_ctrl/OUT/19350101.ocean_month.nc', decode_times = False)
ds1 = ds.zos
ds1 = ds1.to_dataset(name='zos')

dsgrid = xr.open_dataset('/archive/milicak/MOM6-examples/Projects/patchy_NA/work_ctrl/ocean_geometry.nc',
               decode_times = False)

dsgrid = dsgrid.rename_dims({'lonq': 'xq', 'latq': 'yq', 'lonh': 'xh', 'lath': 'yh'})

ds = xr.merge([ds1, dsgrid])
ds = ds.drop_vars({'lonh', 'lonq', 'lath', 'latq'})

ds_arc = xr.open_dataset('/archive/milicak/dataset/MOM6/Arctic_GFDL/ocean_geometry.nc')

# cut the domain for the Arctic
ds_cut = ds.sel(xh=slice(-300, 60), yh=slice(35,91),
                xq=slice(-300, 60), yq=slice(35,91))

# plot it to see
ds_cut['zos'].isel(time=0).plot(figsize=[8, 6], vmin=-2, vmax=2, cmap='needJet2')

# convert lon between -180 to 180
ds_cut = convert_lon(ds_cut)

regrid_domain = xe.Regridder(ds_cut.rename({'geolon': 'lon', 'geolat': 'lat'}),
                             ds_arc.rename({'geolon': 'lon', 'geolat': 'lat'}), 'bilinear',
                                periodic=False, filename='regrid_domain.nc',
                            reuse_weights=True)


# extrapolate ocean values into the land
drowned_ssh = flood_kara(ds_cut['zos'], xdim='xh', ydim='yh')
# plot the new dataset
drowned_ssh.isel(time=0, z=0).plot(figsize=[8, 6], vmin=-2, vmax=2, cmap='jet')

# interpoalte to the regional domain
ssh_ic_arc = regrid_domain(drowned_ssh)

# plot the surface interpolated data
m = Basemap(projection='npstere',boundinglat=40,lon_0=0,resolution='l')
m.fillcontinents(color='grey');
m.pcolormesh(np.copy(ds_arc.geolon),np.copy(ds_arc.geolat),ssh_ic_arc[0,0,:,:],
             latlon=True);plt.colorbar();


# get the last year of the temperature field
ssh_IC_arc = ssh_ic_arc[-12,0,:,:]
ssh_IC_arc = ssh_IC_arc.drop('z')
# create xarray dataset
sshnew_remapped = xr.Dataset()
sshnew_remapped['SSH'] = xr.DataArray(data=ssh_IC_arc.data, dims=('yh', 'xh'))

# one more time to extrapolate to remove nans
sshnew = flood_kara(sshnew_remapped['SSH'], xdim='xh', ydim='yh')
sshnew = sshnew.isel(time=0)
sshnew = sshnew.drop('z')
sshnew = sshnew[0,:,:]
sshnew = sshnew.to_dataset(name='SSH')
sshnew['yh'] = np.arange(0,sshnew.SSH.shape[0])*1.0
sshnew['xh'] = np.arange(0,sshnew.SSH.shape[1])*1.0

dd = sshnew.where(sshnew!=0,np.nan)
# one more time to extrapolate to remove nans
sshnew = flood_kara(dd['SSH'], xdim='xh', ydim='yh')
sshnew = sshnew.isel(time=0)
sshnew = sshnew.drop('z')
sshnew = sshnew[0,:,:]
sshnew = sshnew.to_dataset(name='SSH')
sshnew['yh'] = np.arange(0,sshnew.SSH.shape[0])*1.0
sshnew['xh'] = np.arange(0,sshnew.SSH.shape[1])*1.0
sshnew=sshnew.drop('time')

sshnew.to_netcdf('/archive/milicak/dataset/MOM6/Arctic_GFDL/eta_IC.nc')

