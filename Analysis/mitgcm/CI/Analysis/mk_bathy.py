import os
import numpy as np
import xarray as xr
import scipy.io
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.colors as colors
from scipy.signal import medfilt2d
import netCDF4
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.path import Path
#for interpolation
from scipy.spatial import cKDTree


def lon_lat_to_cartesian(lon, lat, R = 1):
    """
    calculates lon, lat coordinates of a point on a sphere with
    radius R
    """
    lon_r = np.radians(lon)
    lat_r = np.radians(lat)

    x =  R * np.cos(lat_r) * np.cos(lon_r)
    y = R * np.cos(lat_r) * np.sin(lon_r)
    z = R * np.sin(lat_r)
    return x,y,z





inputname = '/media/milicak/DATA1/datasets/world_bathy/bathy_marmara.nc'
# inputname = 'uTSS_lobc_chunk_0365.nos.nc'
outputfile = 'Marmara_bathy_Nx_900_Ny_600.bin'
romsgrdname = '/home/milicak/Analysis/mitgcm/CI/Analysis/ROMS_new_coastline.nc'

ds = xr.open_dataset(inputname)

depth_input = np.copy(ds.z)
lon_input = np.copy(ds.x)
lat_input = np.copy(ds.y)
depth_input = depth_input[::2,::2]
lon_input = lon_input[::2]
lat_input = lat_input[::2]
lon_input, lat_input = np.meshgrid(lon_input,lat_input)

# there is an offset in Emin's lon,lat,bathy files and this trickled down to
# Sannino's and also Shyfem models
# depth_input = np.copy(ds.total_depth)
# lon_input = np.copy(ds.longitude)
# lat_input = np.copy(ds.latitude)

df = xr.open_dataset(romsgrdname)
lon_roms = np.copy(df.lon_rho)
lat_roms = np.copy(df.lat_rho)

nj,ni = df.lon_rho.shape

# 1st Way
aa = np.column_stack((lon_input.flatten(),lat_input.flatten()))
tree = cKDTree(aa)
d, inds = tree.query(np.transpose(np.vstack((np.copy(df.lon_rho).flatten(),
                                             np.copy(df.lat_rho).flatten()))),
                     k = 2)
w = 1.0 / d**2
depth_mitgcm = np.sum(w * depth_input.flatten()[inds], axis=1) / np.sum(w, axis=1)
depth_mitgcm.shape = df.lon_rho.shape

dnm = np.copy(depth_mitgcm*df.mask_rho)
dnm[dnm<0] = 0
dnm[np.logical_and(dnm>0,dnm<4)]=4
depth_mitgcm = dnm
depth_mitgcm.astype('>f8').tofile(outputfile)



# 2nd Way
xs, ys, zs = lon_lat_to_cartesian(lon_input.flatten(), lat_input.flatten())
del zs
xt, yt, zt = lon_lat_to_cartesian(lon_roms.flatten(), lat_roms.flatten())
del zt
bb = np.column_stack((xs,ys))
tree2 = cKDTree(bb)
d2, inds2 = tree2.query(np.transpose(np.vstack((xt,yt))),k = 2)
w2 = 1.0 / d2**2
depth2_mitgcm = np.sum(w2 * depth_input.flatten()[inds2], axis=1) / np.sum(w2, axis=1)
depth2_mitgcm.shape = df.lon_rho.shape

# 3rd Way
# m = Basemap(llcrnrlon=26.0,llcrnrlat=40.0,urcrnrlon=30.,urcrnrlat=41.4,
#             rsphere=(6378137.00,6356752.3142),
#             resolution='f',projection='merc',
#             lat_0=40.,lon_0=20.,lat_ts=20.)
# m.drawcoastlines(linewidth=0.2)
#
# longitude,latitude = m(np.copy(ds.longitude)-0.004,np.copy(ds.latitude)-0.0015)
# longitude2,latitude2 = m(np.copy(df.lon_rho),np.copy(df.lat_rho))
# aa = np.column_stack((longitude.flatten(),latitude.flatten()))
# tree = cKDTree(aa)
# d, inds = tree.query(np.transpose(np.vstack((np.copy(longitude2).flatten(),
#                                              np.copy(latitude2).flatten()))),
#                      k = 4)
# w = 1.0 / d**2
# depth_mitgcm = np.sum(w * depth_shyfem.flatten()[inds], axis=1) / np.sum(w, axis=1)
# depth_mitgcm.shape = df.lon_rho.shape
# # plt.tripcolor(longitude,latitude,ds.element_index-1,ds.total_depth,cmap='needJet2',shading='gouraud');plt.colorbar();
