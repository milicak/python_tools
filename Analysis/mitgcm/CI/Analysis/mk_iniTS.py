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



inputname = 'uTSS_lobc_chunk_0365.nos.nc'
outputfilet = 'Tinit_Nx_900_Ny_600.bin'
outputfiles = 'Sinit_Nx_900_Ny_600.bin'
romsgrdname = '/home/milicak/Analysis/mitgcm/CI/Analysis/ROMS_new_coastline.nc'
mitgrdname = 'Marmara_bathy_Nx_900_Ny_600.nc'

gr = xr.open_dataset(mitgrdname)

ds = xr.open_dataset(inputname)

temp_input = np.copy(ds.temperature[0,:,:])
salt_input = np.copy(ds.salinity[0,:,:])
lon_input = np.copy(ds.longitude)
lat_input = np.copy(ds.latitude)


df = xr.open_dataset(romsgrdname)
lon_roms = np.copy(df.lon_rho)
lat_roms = np.copy(df.lat_rho)

nj,ni = df.lon_rho.shape

aa = np.column_stack((lon_input.flatten(),lat_input.flatten()))
tree = cKDTree(aa)
d, inds = tree.query(np.transpose(np.vstack((np.copy(df.lon_rho).flatten(),
                                             np.copy(df.lat_rho).flatten()))),
                     k = 2)
w = 1.0 / d**2

temp_mitgcm = np.zeros((93,nj,ni))
salt_mitgcm = np.zeros((93,nj,ni))
for kind in range(0,93):
    kind
    tmp = temp_input[:,kind]
    tmp_mitgcm = np.sum(w * tmp.flatten()[inds], axis=1) / np.sum(w, axis=1)
    tmp_mitgcm.shape = df.lon_rho.shape
    temp_mitgcm[kind,:,:] = tmp_mitgcm
    tmp = salt_input[:,kind]
    tmp_mitgcm = np.sum(w * tmp.flatten()[inds], axis=1) / np.sum(w, axis=1)
    tmp_mitgcm.shape = df.lon_rho.shape
    salt_mitgcm[kind,:,:] = tmp_mitgcm





# Write to a file size should be Nz, Ny, Nx
temp_mitgcm.astype('>f8').tofile(outputfilet)
salt_mitgcm.astype('>f8').tofile(outputfiles)





