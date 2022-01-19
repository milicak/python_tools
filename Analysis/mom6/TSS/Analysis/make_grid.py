# This file is designed to be cut and pasted into an ipython --pylab
# session. Otherwise, you'll need to "import np as np" then
# convert "array" to "np.array".
import os
import numpy as np
import scipy.io
from mpl_toolkits.basemap import Basemap, shiftgrid
import matplotlib.colors as colors
from scipy.signal import medfilt2d
import netCDF4
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import pyroms
import pyroms_toolbox
from bathy_smoother import *
from matplotlib.path import Path
#for interpolation
from scipy.spatial import cKDTree


def angle_p1p2(p1, p2):
    """Angle at center of sphere between two points on the surface of the sphere.
    Positions are given as (latitude,longitude) tuples measured in degrees."""
    phi1 = np.deg2rad( p1[0] )
    phi2 = np.deg2rad( p2[0] )
    dphi_2 = 0.5 * ( phi2 - phi1 )
    dlambda_2 = 0.5 * np.deg2rad( p2[1] - p1[1] )
    a = np.sin( dphi_2 )**2 + np.cos( phi1 ) * np.cos( phi2 ) * ( np.sin( dlambda_2 )**2 )
    c = 2. * np.arctan2( np.sqrt(a), np.sqrt( 1. - a ) )
    return c

def spherical_angle(v1, v2, v3):
    """Returns angle v2-v1-v3 i.e betweeen v1-v2 and v1-v3."""
    # vector product between v1 and v2
    px = v1[1]*v2[2] - v1[2]*v2[1]
    py = v1[2]*v2[0] - v1[0]*v2[2]
    pz = v1[0]*v2[1] - v1[1]*v2[0]
    # vector product between v1 and v3
    qx = v1[1]*v3[2] - v1[2]*v3[1]
    qy = v1[2]*v3[0] - v1[0]*v3[2]
    qz = v1[0]*v3[1] - v1[1]*v3[0]

    ddd = (px*px+py*py+pz*pz)*(qx*qx+qy*qy+qz*qz)
    ddd = (px*qx+py*qy+pz*qz) / np.sqrt(ddd)
    angle = np.arccos( ddd );
    return angle
def spherical_quad(lat,lon):
    """Returns area of spherical quad (bounded by great arcs)."""
    # x,y,z are 3D coordinates
    d2r = np.deg2rad(1.)
    x = np.cos(d2r*lat)*np.cos(d2r*lon)
    y = np.cos(d2r*lat)*np.sin(d2r*lon)
    z = np.sin(d2r*lat)
    c0 = (x[:-1,:-1],y[:-1,:-1],z[:-1,:-1])
    c1 = (x[:-1,1:],y[:-1,1:],z[:-1,1:])
    c2 = (x[1:,1:],y[1:,1:],z[1:,1:])
    c3 = (x[1:,:-1],y[1:,:-1],z[1:,:-1])
    a0 = spherical_angle( c1, c0, c2)
    a1 = spherical_angle( c2, c1, c3)
    a2 = spherical_angle( c3, c2, c0)
    a3 = spherical_angle( c0, c3, c1)
    return a0+a1+a2+a3-2.*np.pi



# generate the mask
#for verts in map.coastsegs:
#    hgrd.mask_polygon(verts)
# alternate version from johan.navarro.padron

hgrd = scipy.io.loadmat('MITgcm_TSS_grid.mat')
lon_rho = np.transpose(hgrd['lon_rho']) #[1:-1,1:-1])
lon_psi = np.transpose(hgrd['lon_psi']) #[1:-2,1:-2])
lat_rho = np.transpose(hgrd['lat_rho']) #[1:-1,1:-1])
lat_psi = np.transpose(hgrd['lat_psi']) #[1:-2,1:-2])

lon_u = np.transpose(hgrd['lon_u'])
lon_v = np.transpose(hgrd['lon_v'])
lat_u = np.transpose(hgrd['lat_u'])
lat_v = np.transpose(hgrd['lat_v'])

nj,ni = lon_rho.shape
nj -=2; ni -=2
print('nj=%i, nj=%i'%(nj,ni))

# sni - Number of columns in x-direction
# snj - Number of rows in y-direction
# x[snj+1,sni+1] - x-coordinate of the grid vertices
# y[snj+1,sni+1] - y-coordinate of the grid vertices
# area[snj,sni] - area of each cell
# dx[snj+1,sni] - length of theedges of constant j-index
# dy[snj,sni+1] - length of the edges of constant -index

# Supergrid shape
snj,sni = 2*np.array([nj,ni]) # Smallest useful super-grid has a multiplier of 2

# Declare shapes
lon = np.zeros((snj+1,sni+1))
lat = np.zeros((snj+1,sni+1))
area = np.zeros((snj,sni))
dx = np.zeros((snj+1,sni))
dy = np.zeros((snj,sni+1))
angle = np.zeros((snj+1,sni+1))

# Copy in data from ROMS file
lon[::2,::2] = lon_psi[:,:] # Cell corners
lon[1::2,1::2] = lon_rho[1:-1,1:-1] # Cell centers (drop outside row and column)
lon[1::2,::2] = lon_u[1:-1,:] # U-points (drop outside row)
lon[::2,1::2] = lon_v[:,1:-1] # V-points (drop outside column)
lat[::2,::2] = lat_psi[:,:] # Cell corners
lat[1::2,1::2] =lat_rho[1:-1,1:-1] # Cell centers (drop outside row and column)
lat[1::2,::2] = lat_u[1:-1,:] # U-points (drop outside row)
lat[::2,1::2] = lat_v[:,1:-1] # V-points (drop outside column)

# tApproximate edge lengths as great arcs
R = 6370.e3 # Radius of sphere
dx[:,:] = R*angle_p1p2( (lat[:,1:],lon[:,1:]), (lat[:,:-1],lon[:,:-1]) )
dy[:,:] = R*angle_p1p2( (lat[1:,:],lon[1:,:]), (lat[:-1,:],lon[:-1,:]) )

# Approximate angles using centered differences in interior
angle[:,1:-1] = np.arctan( (lat[:,2:]-lat[:,:-2]) /
                          ((lon[:,2:]-lon[:,:-2])*np.cos(np.deg2rad(lat[:,1:-1]))) )
# Approximate angles using side differences on left/right edges
angle[:,0] = np.arctan( (lat[:,1]-lat[:,0]) / ((lon[:,1]-lon[:,0])*np.cos(np.deg2rad(lat[:,0]))) )
angle[:,-1] = np.arctan( (lat[:,-1]-lat[:,-2]) /
                        ((lon[:,-1]-lon[:,-2])*np.cos(np.deg2rad(lat[:,-1]))) )

area = dx[:-1,:]*dy[:,:-1]
# area = area[:-1,:-1]

# Create a mosaic file
rg = scipy.io.netcdf_file('/okyanus/users/milicak/dataset/MOM6/TSS/ocean_hgrid.nc','w')
# Dimensions
rg.createDimension('nx',sni)
rg.createDimension('nxp1',sni+1)
rg.createDimension('ny',snj)
rg.createDimension('nyp1',snj+1)
# rg.createDimension('string',255)
rg.createDimension('string',5)
# Variables
hx = rg.createVariable('x','float32',('nyp1','nxp1',))
hx.units = 'degrees'
hy = rg.createVariable('y','float32',('nyp1','nxp1',))
hy.units = 'degrees'
hdx = rg.createVariable('dx','float32',('nyp1','nx',))
hdx.units = 'meters'
hdy = rg.createVariable('dy','float32',('ny','nxp1',))
hdy.units = 'meters'
harea = rg.createVariable('area','float32',('ny','nx',))
harea.units = 'meters^2'
hangle = rg.createVariable('angle_dx','float32',('nyp1','nxp1',))
hangle.units = 'degrees'
htile = rg.createVariable('tile','c',('string',))
# Values
hx[:] = lon
hy[:] = lat
hdx[:] = dx
hdy[:] = dy
harea[:] = area
hangle[:] = angle
htile[:5] = 'tile1'
rg.close()
