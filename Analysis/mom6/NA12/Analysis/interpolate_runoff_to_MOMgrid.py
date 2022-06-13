import numpy as np
import xesmf
import xarray

root_folder = '~/dataset/MOM6/NA12/'
out_file = root_folder + 'glofas-era5_NA12_gridded_1996.nc'

def get_coast_mask(mask_file):
    mask = xarray.open_dataset(mask_file)

    # Alistair's method of finding coastal cells
    ocn_mask = mask['mask'].values
    cst_mask = 0 * ocn_mask # All land should be 0
    is_ocean = ocn_mask > 0
    cst_mask[(is_ocean) & (np.roll(ocn_mask, 1, axis=1) == 0)] = 1 # Land to the west
    cst_mask[(is_ocean) & (np.roll(ocn_mask, -1, axis=1) == 0)] = 1 # Land to the east
    cst_mask[(is_ocean) & (np.roll(ocn_mask, 1, axis=0) == 0)] = 1 # Land to the south
    cst_mask[(is_ocean) & (np.roll(ocn_mask, -1, axis=0) == 0)] = 1 # Land to the north

    # Model boundaries are not coasts
    cst_mask[0, :] = 0
    cst_mask[:, 0] = 0
    cst_mask[-1, :] = 0
    cst_mask[:, -1] = 0

    return cst_mask


glofas = xr.open_dataset('~/dataset/MOM6/NA12/glofas-era5_NA12_1996.nc')
hgrid = xr.open_dataset('~/dataset/MOM6/NA12/ocean_hgrid.nc')
coast_mask = get_coast_mask('/okyanus/users/milicak/dataset/MOM6/NA12/ocean_mask.nc')

glofas_latb = np.arange(glofas['lat'][0]+.05, glofas['lat'][-1]-.051, -.1)
glofas_lonb = np.arange(glofas['lon'][0]-.05, glofas['lon'][-1]+.051, .1)
lon = hgrid.x[1::2, 1::2]
lonb = hgrid.x[::2, ::2]
lat = hgrid.y[1::2, 1::2]
latb = hgrid.y[::2, ::2]
# From Alistair
area = (hgrid.area[::2, ::2] + hgrid.area[1::2, 1::2]) + (hgrid.area[1::2, ::2] + hgrid.area[::2, 1::2])

# Conservatively interpolate runoff onto MOM grid
glofas_to_mom_con = xesmf.Regridder(
    {'lon': glofas.lon, 'lat': glofas.lat, 'lon_b': glofas_lonb, 'lat_b': glofas_latb},
    {'lat': lat, 'lon': lon, 'lat_b': latb, 'lon_b': lonb},
    method='conservative',
    periodic=True,
    reuse_weights=True
)

# Interpolate only from GloFAS points that are river end points.
# glofas_regridded = glofas_to_mom_con(glofas)
glofas_regridded = glofas_to_mom_con(glofas.runoff.fillna(0))
glofas_regridded = glofas_regridded.rename({'nyp': 'ny', 'nxp': 'nx'}).values


# Flatten mask and coordinates to 1D
flat_mask = coast_mask.ravel().astype('bool')
coast_lon = lon.values.ravel()[flat_mask]
coast_lat = lat.values.ravel()[flat_mask]
mom_id = np.arange(np.prod(coast_mask.shape))

# Use xesmf to find the index of the nearest coastal cell
# for every grid cell in the MOM domain
coast_to_mom = xesmf.Regridder(
    {'lat': coast_lat, 'lon': coast_lon},
    {'lat': lat, 'lon': lon},
    method='nearest_s2d',
    locstream_in=True,
    reuse_weights=True
)

coast_id = mom_id[flat_mask]
nearest_coast = coast_to_mom(coast_id).ravel()

# Raw runoff on MOM grid, reshaped to 2D (time, grid_id)
raw = glofas_regridded.reshape([glofas_regridded.shape[0], -1])

 # Zero array that will be filled with runoff at coastal cells
filled = np.zeros_like(raw)

# Loop over each coastal cell and fill the result array
# with the sum of runoff for every grid cell that
# has this coastal cell as its closest coastal cell
for i in coast_id:
    filled[:, i] = raw[:, nearest_coast == i].sum(axis=1)

# Reshape back to 3D
filled_reshape = filled.reshape(glofas_regridded.shape)

# Convert to xarray
ds = xarray.Dataset({
    'runoff': (['time', 'y', 'x'], filled_reshape),
    'area': (['y', 'x'], area.data),
    'lat': (['y', 'x'], lat.data),
    'lon': (['y', 'x'], lon.data)
    },
    coords={'time': glofas['time'].data, 'y': np.arange(filled_reshape.shape[1]), 'x': np.arange(filled_reshape.shape[2])}
)

# Drop '_FillValue' from all variables when writing out
all_vars = list(ds.data_vars.keys()) + list(ds.coords.keys())
encodings = {v: {'_FillValue': None} for v in all_vars}

# Make sure time has the right units and datatype
# otherwise it will become an int and MOM will fail.
encodings['time'].update({
    'units': 'days since 1950-01-01',
    'dtype': np.float,
    'calendar': 'gregorian'
})

ds['time'].attrs = {'cartesian_axis': 'T'}
ds['x'].attrs = {'cartesian_axis': 'X'}
ds['y'].attrs = {'cartesian_axis': 'Y'}
ds['lat'].attrs = {'units': 'degrees_north'}
ds['lon'].attrs = {'units': 'degrees_east'}
ds['runoff'].attrs = {'units': 'kg m-2 s-1'}

# Write out
ds.to_netcdf(
    out_file,
    unlimited_dims=['time'],
    format='NETCDF3_64BIT',
    encoding=encodings,
    engine='netcdf4'
)
ds.close()


