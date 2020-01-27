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




inputname = 'ssh.nc'
outputfile = 'SSH_INI_BS_Turkey_coast.bin'
romsgrdname = 'ROMS_BS_forecast_Turkey_coast.nc'

ds = xr.open_dataset(inputname);
ssh_nemo = np.copy(ds.sossheig[0,:,:])
lon_nemo = np.copy(ds.lon);
lat_nemo = np.copy(ds.lat);
lon_nemo,lat_nemo = np.meshgrid(lon_nemo,lat_nemo)

df = xr.open_dataset(romsgrdname)
lon_roms = np.copy(df.lon_rho)
lat_roms = np.copy(df.lat_rho)

nj,ni = df.lon_rho.shape

# Use cKDTree
# xs, ys, zs = lon_lat_to_cartesian(lon_gebco.flatten(), lat_gebco.flatten())
# del zs
# xt, yt, zt = lon_lat_to_cartesian(lon_roms.flatten(), lat_roms.flatten())
# del zt
# bb = np.column_stack((xs,ys))
aa = np.column_stack((lon_nemo.flatten(),lat_nemo.flatten()))
tree = cKDTree(aa)

d, inds = tree.query(np.transpose(np.vstack((np.copy(df.lon_rho).flatten(),
                                             np.copy(df.lat_rho).flatten()))),
                     k = 4)
w = 1.0 / d**2

ssh_mitgcm = np.sum(w * ssh_nemo.flatten()[inds], axis=1) / np.sum(w, axis=1)
ssh_mitgcm.shape = df.lon_rho.shape
mask = np.isnan(ssh_mitgcm)
ssh_mitgcm[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), ssh_mitgcm[~mask])

# Write to a file
ssh_mitgcm.astype('>f4').tofile(outputfile)

# aa=np.zeros((80,84))
# aa[np.array(np.isnan(ssh_mitgcm), dtype=bool)==1 & (np.copy(df.mask_rho)==1)]=1
