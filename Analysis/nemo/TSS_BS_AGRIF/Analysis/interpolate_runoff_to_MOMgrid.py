import numpy as np
from matplotlib.path import Path as mpPath
from mpl_toolkits.basemap import Basemap
import numpy.ma as ma
import glob
import xarray as xr
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import xesmf
import xarray
import os

root_folder = '~/dataset/NEMO/TSS_BS_AGRIF/'

def get_coast_mask(mask_file):
    df = xr.open_dataset(mask_file)
    mask = xr.where(df.bathy_metry[0,:,:]!=0,1,0)

    # Alistair's method of finding coastal cells
    ocn_mask = mask.values
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


hgrid = xr.open_dataset('~/dataset/NEMO/TSS_BS_AGRIF/domain_cfg.nc')
coast_mask = get_coast_mask('~/dataset/NEMO/TSS_BS_AGRIF/domain_cfg.nc')

year1 = 1996
fname = '~/dataset/GLOFAS/glofas-era5_' + str(year1) + '.nc'
df = xr.open_dataset(fname)
df['lon']=df.lon-360
df = df.reindex(lat=list(reversed(df.lat)))
# Convert m3/s to kg/m2/s
# Borrowed from https://xgcm.readthedocs.io/en/latest/xgcm-examples/05_autogenerate.html
distance_1deg_equator = 111000.0
dlon = dlat = 0.1  # GloFAS grid spacing
dx = dlon * np.cos(np.deg2rad(df.lat)) * distance_1deg_equator
dy = ((df.lon * 0) + 1) * dlat * distance_1deg_equator
glofas_area = dx * dy

uparea = xr.open_dataarray('/okyanus/users/milicak/dataset/MOM6/NA12/upArea.nc')
# Find river end points by looking for local maxima in upstream area.
uparea = uparea.fillna(0).values
points = np.zeros_like(uparea)
window = 2  # look with +- this number of grid points
ni, nj = uparea.shape
for i in range(window, ni-window):
    for j in range(window, nj-window):
        sub = uparea[i-window:i+window+1, j-window:j+window+1]
        point = uparea[i, j]
        # A river end point has a reasonably large upstream area
        # and is a local maximum
        if point > 1e6 and sub.max() == point:
            points[i, j] = 1

points = np.flipud(points)
lon = np.copy(df.lon)
lat = np.copy(df.lat)
lon_mom = np.copy(hgrid.nav_lon)
lat_mom = np.copy(hgrid.nav_lat)

index = np.abs(lon-lon_mom.min()).argmin()
ind1 = index-1
index = np.abs(lon-lon_mom.max()).argmin()
ind2 = index+1
index = np.abs(lat-lat_mom.min()).argmin()
ind3 = index
index = np.abs(lat-lat_mom.max()).argmin()
ind4 = index+1
inside = np.zeros((df.dis24.shape[1],df.dis24.shape[2]))
inside[ind3:ind4+1,ind1:ind2+1] = 1

# if os.path.exists('inside_matrix.mat'):
#     print('mehmet')
# else:
#     m = Basemap(projection='cea',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='c')
#     lon = np.copy(df.lon)
#     lat = np.copy(df.lat[::-1])
#     lon, lat = np.meshgrid(lon,lat)
#     lon1, lat1 = m(lon,lat)
#     lon_mom = np.copy(hgrid.nav_lon)
#     lat_mom = np.copy(hgrid.nav_lat)
#     lon_mom1, lat_mom1 = m(lon_mom,lat_mom)
#     vertices = np.transpose(np.array([lon_mom1.flatten(),lat_mom1.flatten()]))
#     pnts = np.transpose(np.array([lon1.flatten(),lat1.flatten()]))
#     path = mpPath(vertices)
#     inside = path.contains_points(pnts)
#     ny = lon.shape[0]
#     nx = lon.shape[1]
#     inside = np.reshape(inside,((ny,nx)))
#     inside = np.double(inside)
#
lon = np.copy(hgrid.glamt[0,:,:])
lat = np.copy(hgrid.gphit[0,:,:])
lonb = np.zeros((lon.shape[0]+1,lon.shape[1]+1))
latb = np.zeros((lon.shape[0]+1,lon.shape[1]+1))
lonb[1:,1:] = np.copy(hgrid.glamf[0,:,:])
latb[1:,1:] = np.copy(hgrid.gphif[0,:,:])
# 0.02 degrees
lonb[1:,0] = lonb[1:,1]-0.02
latb[0,1:] = latb[1,1:]-0.02
lonb[0,:] = lonb[1,:]
latb[:,0] = latb[:,1]
area = np.copy(hgrid.e1t[0,:,:]*hgrid.e2t[0,:,:])

for year in range(1993,1995):
    print(year)
    fname = '~/dataset/GLOFAS/glofas-era5_' + str(year) + '.nc'
    out_file = root_folder + 'glofas-era5_TSS_BS_AGRIF_' + str(year) + '.nc'
    glofas_tmp = xr.open_dataset(fname)
    glofas_tmp['lon'] = glofas_tmp.lon-360
    glofas_tmp = glofas_tmp.reindex(lat=list(reversed(glofas_tmp.lat)))
    # convert m3/s to kgm-2s-1
    glofas = glofas_tmp * (1000.0*inside*points) / glofas_area
    print(glofas)
    # glofas = glofas*inside*points

    glofas_latb = np.arange(glofas['lat'][0]-.05, glofas['lat'][-1]+.051, .1)
    glofas_lonb = np.arange(glofas['lon'][0]-.05, glofas['lon'][-1]+.051, .1)

    # Conservatively interpolate runoff onto MOM grid
    # print({'lon': glofas.lon, 'lat': glofas.lat, 'lon_b': glofas_lonb, 'lat_b': glofas_latb})
    print(lon.shape)
    print(lonb.shape)
    print(lat.shape)
    print(latb.shape)
    glofas_to_nemo_con = xesmf.Regridder(
        {'lon': glofas.lon, 'lat': glofas.lat, 'lon_b': glofas_lonb, 'lat_b': glofas_latb},
        {'lat': lat, 'lon': lon, 'lat_b': latb, 'lon_b': lonb},
        method='conservative',
        periodic=True,
        reuse_weights=True
    )

    # Interpolate only from GloFAS points that are river end points.
    # glofas_regridded = glofas_to_nemo_con(glofas)
    glofas_regridded = glofas_to_nemo_con(glofas.dis24.fillna(0))
    glofas_regridded = np.copy(glofas_regridded)
    # glofas_regridded = glofas_regridded.rename({'nyp': 'ny', 'nxp': 'nx'}).values

    # Flatten mask and coordinates to 1D
    flat_mask = coast_mask.ravel().astype('bool')
    coast_lon = lon.ravel()[flat_mask]
    coast_lat = lat.ravel()[flat_mask]
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
        'runoff': (['time_counter', 'y', 'x'], filled_reshape),
        'area': (['y', 'x'], area.data),
        'lat': (['y', 'x'], lat.data),
        'lon': (['y', 'x'], lon.data)
        },
        coords={'time_counter': glofas['time'].data, 'y': np.arange(filled_reshape.shape[1]), 'x': np.arange(filled_reshape.shape[2])}
    )

    # Drop '_FillValue' from all variables when writing out
    all_vars = list(ds.data_vars.keys()) + list(ds.coords.keys())
    encodings = {v: {'_FillValue': None} for v in all_vars}
    # Make sure time has the right units and datatype
    # otherwise it will become an int and MOM will fail.
    encodings['time_counter'].update({
        'units': 'days since 1950-01-01',
        'dtype': np.float,
        'calendar': 'gregorian'
    })

    ds['time_counter'].attrs = {'cartesian_axis': 'T'}
    ds['x'].attrs = {'cartesian_axis': 'X'}
    ds['y'].attrs = {'cartesian_axis': 'Y'}
    ds['lat'].attrs = {'units': 'degrees_north'}
    ds['lon'].attrs = {'units': 'degrees_east'}
    ds['runoff'].attrs = {'units': 'kg m-2 s-1'}

    # Write out
    ds.to_netcdf(
        out_file,
        unlimited_dims=['time_counter'],
        format='NETCDF3_64BIT',
        encoding=encodings,
        engine='netcdf4'
    )
    ds.close()




