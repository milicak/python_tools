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




inputname = 'dnm_temp.nc'
outputfile = 'TEMP_INI_BS_Turkey_coast.bin'
romsgrdname = 'ROMS_BS_forecast_Turkey_coast.nc'

ds = xr.open_dataset(inputname);
var_nemo = np.copy(ds.votemper[0,:,:,:])
lon_nemo = np.copy(ds.lon);
lat_nemo = np.copy(ds.lat);
lon_nemo,lat_nemo = np.meshgrid(lon_nemo,lat_nemo)
zr_nemo = np.concatenate(([0],ds.depth))

df = xr.open_dataset(romsgrdname)
lon_roms = np.copy(df.lon_rho)
lat_roms = np.copy(df.lat_rho)

nj,ni = df.lon_rho.shape

# Use cKDTree
aa = np.column_stack((lon_nemo.flatten(),lat_nemo.flatten()))
tree = cKDTree(aa)

d, inds = tree.query(np.transpose(np.vstack((np.copy(df.lon_rho).flatten(),
                                             np.copy(df.lat_rho).flatten()))),
                     k = 4)
w = 1.0 / d**2


var3d = np.zeros((var_nemo.shape[0]+1,nj,ni))
for kind in range(0,var_nemo.shape[0]):
    var_mitgcm = np.sum(w * var_nemo[kind,:,:].flatten()[inds], axis=1) / np.sum(w, axis=1)
    var_mitgcm.shape = df.lon_rho.shape
    mask = np.isnan(var_mitgcm)
    if (ma.all(mask)):
        var_mitgcm[mask] = 5.0;
    else:
        var_mitgcm[mask] = np.interp(np.flatnonzero(mask),
                                     np.flatnonzero(~mask), var_mitgcm[~mask])


    if(kind==0):
        var3d[kind,:,:] = var_mitgcm
        var3d[kind+1,:,:] = var_mitgcm
    else:
        var3d[kind+1,:,:] = var_mitgcm


# for iind in range(0,ni):
#     for jind in range(0,nj):



# Linear interpolation for the MITgcm vertical grid

# Write to a file size should be Nz, Ny, Nx
# ssh_mitgcm.astype('>f4').tofile(outputfile)

