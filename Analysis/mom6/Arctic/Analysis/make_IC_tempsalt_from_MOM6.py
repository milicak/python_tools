from vcr import utils, conserve
from HCtFlood.kara import flood_kara
import xesmf as xe
# import GFDL_xr


def convert_lon(ds, lon='xh'):
    ds.coords[lon] = (ds.coords[lon] + 180) % 360 - 180
    ds = ds.sortby(ds[lon])
    return ds


ds = xr.open_dataset('/archive/milicak/MOM6-examples/Projects/patchy_NA/work_ctrl/OUT/19350101.ocean_month_z.nc', decode_times = False)
# ds = GFDL_xr.open_dataset('/archive/milicak/MOM6-examples/Projects/patchy_NA/work_ctrl/OUT/19350101.ocean_month_z.nc', decode_times = False)
ds = ds.drop_vars({'uo', 'vo', 'umo', 'vmo', 'obvfsq', 'agessc', 'uhml', 'vhml',
              'uhGM', 'vhGM', 'uh', 'vh', 'T_adx', 'T_ady', 'S_adx', 'S_ady',
              'difvho', 'difvso', 'Kd_interface', 'Kd_shear', 'Kd_itides',
              'Kd_BBL', 'Kd_ePBL', 'average_T1', 'average_T2', 'average_DT'})

dsgrid = xr.open_dataset('/archive/milicak/MOM6-examples/Projects/patchy_NA/work_ctrl/ocean_geometry.nc',
               decode_times = False)

dsgrid = dsgrid.rename_dims({'lonq': 'xq', 'latq': 'yq', 'lonh': 'xh', 'lath': 'yh'})

ds = xr.merge([ds, dsgrid])
ds = ds.drop_vars({'lonh', 'lonq', 'lath', 'latq'})

ds_arc = xr.open_dataset('/archive/milicak/dataset/MOM6/Arctic_GFDL/ocean_geometry.nc')

# vertical grid for the regional model
vsgrid = xr.open_dataset('/archive/milicak/dataset/MOM6/Arctic_GFDL/vgrid_75_2m.nc')
zw = np.zeros(len(vsgrid.dz) + 1)
zw[1:] = np.cumsum(vsgrid.dz)
depth_tgt = 0.5*(zw[0:-1]+zw[1::])
depth_bnds_tgt = zw

# cut the domain for the Arctic
ds_cut = ds.sel(xh=slice(-300, 60), yh=slice(35,91),
                xq=slice(-300, 60), yq=slice(35,91))

# plot it to see
ds_cut['thetao'].isel(time=0, z_l=0).plot(figsize=[8, 6], vmin=-2, vmax=20, cmap='jet')

# convert lon between -180 to 180
ds_cut = convert_lon(ds_cut)

regrid_domain = xe.Regridder(ds_cut.rename({'geolon': 'lon', 'geolat': 'lat'}),
                             ds_arc.rename({'geolon': 'lon', 'geolat': 'lat'}), 'bilinear',
                                periodic=False, filename='regrid_domain.nc',
                            reuse_weights=True)


# extrapolate ocean values into the land
drowned_temp = flood_kara(ds_cut['thetao'], xdim='xh', ydim='yh', zdim='z_l')
# plot the new dataset
drowned_temp.isel(time=0, z_l=0).plot(figsize=[8, 6], vmin=-2, vmax=20, cmap='jet')

# extrapolate vertically to fill NaN values
alldrowned_temp = drowned_temp.ffill(dim='z_l')

# interpoalte to the regional domain
temp_ic_arc = regrid_domain(alldrowned_temp)

# plot the surface interpolated data
temp_ic_arc.isel(time=0, z_l=0).plot(figsize=[8, 6], x='lon', y='lat', vmin=-2,
                                     vmax=20, cmap='jet')


m = Basemap(projection='npstere',boundinglat=40,lon_0=0,resolution='l')
m.fillcontinents(color='grey');
m.pcolormesh(np.copy(ds_arc.geolon),np.copy(ds_arc.geolat),temp_ic_arc[0,0,:,:],
             latlon=True);plt.colorbar();


# get the last year of the temperature field
temp_IC_arc = temp_ic_arc[-12,:,:]
temp_IC_arc = temp_IC_arc.to_dataset(name='PTEMP')
temp_IC_arc = temp_IC_arc.rename({'time': 'TIME', 'z_l': 'DEPTH', 'lon': 'LON', 'lat': 'LAT'})


# vertical interpolation to the regional grid
# first source vertical grid information
depth_src = np.copy(temp_IC_arc.DEPTH)
depth_bnds_src = np.array([0.000e+00, 5.000e+00, 1.500e+01, 2.500e+01, 4.000e+01, 6.250e+01,
                           8.750e+01, 1.125e+02, 1.375e+02, 1.750e+02, 2.250e+02, 2.750e+02,
                           3.500e+02, 4.500e+02, 5.500e+02, 6.500e+02, 7.500e+02, 8.500e+02,
                           9.500e+02, 1.050e+03, 1.150e+03, 1.250e+03, 1.350e+03, 1.450e+03,
                           1.625e+03, 1.875e+03, 2.250e+03, 2.750e+03, 3.250e+03, 3.750e+03,
                           4.250e+03, 4.750e+03, 5.250e+03, 5.750e+03, 6.250e+03, 6.750e+03])

remapping = conserve.create_remapping_matrix(depth_bnds_src, depth_bnds_tgt, strict=False)
remapping_strict = conserve.create_remapping_matrix(depth_bnds_src, depth_bnds_tgt, strict=True)

# remap to the regional grid
temp_remapped = conserve.vertical_remap_z2z(temp_IC_arc['PTEMP'].values, remapping)
# create xarray dataset
tempnew_remapped = xr.Dataset()
tempnew_remapped['PTEMP'] = xr.DataArray(data=temp_remapped, dims=('z_l', 'yh', 'xh'))
# one more time to extrapolate to remove nans
tempnew = flood_kara(tempnew_remapped['PTEMP'], xdim='xh', ydim='yh', zdim='z_l')
tempnew = tempnew.isel(time=0)
tempnew = tempnew.drop('time')
tempnew = tempnew.to_dataset(name='PTEMP')
tempnew['z_l'] = depth_tgt
tempnew['yh'] = np.arange(0,tempnew.PTEMP.shape[1])*1.0
tempnew['xh'] = np.arange(0,tempnew.PTEMP.shape[2])*1.0

tempnew.to_netcdf('/archive/milicak/dataset/MOM6/Arctic_GFDL/temp_IC.nc')

# extrapolate salinity ocean values into the land
drowned_salt = flood_kara(ds_cut['so'], xdim='xh', ydim='yh', zdim='z_l')
# extrapolate vertically to fill NaN values
alldrowned_salt = drowned_salt.ffill(dim='z_l')
# interpoalte to the regional domain
salt_ic_arc = regrid_domain(alldrowned_salt)
# get the last time step of the temperature field
salt_IC_arc = salt_ic_arc[-12,:,:]
salt_IC_arc = salt_IC_arc.to_dataset(name='SALT')
salt_IC_arc = salt_IC_arc.rename({'time': 'TIME', 'z_l': 'DEPTH', 'lon': 'LON', 'lat': 'LAT'})

salt_remapped = conserve.vertical_remap_z2z(salt_IC_arc['SALT'].values, remapping)
# create xarray dataset
saltnew_remapped = xr.Dataset()
saltnew_remapped['SALT'] = xr.DataArray(data=salt_remapped, dims=('z_l', 'yh', 'xh'))
# one more time to extrapolate to remove nans
saltnew = flood_kara(saltnew_remapped['SALT'], xdim='xh', ydim='yh', zdim='z_l')
saltnew = saltnew.isel(time=0)
saltnew = saltnew.drop('time')
saltnew = saltnew.to_dataset(name='SALT')
saltnew_remapped['z_l'] = depth_tgt
saltnew['yh'] = np.arange(0,saltnew.SALT.shape[1])*1.0
saltnew['xh'] = np.arange(0,saltnew.SALT.shape[2])*1.0

saltnew.to_netcdf('/archive/milicak/dataset/MOM6/Arctic_GFDL/salt_IC.nc')


