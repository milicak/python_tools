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

# Grid dimension
# x-direction
Lp = 11000
# y-direction
Mp = 6060


print('generating the projection')

# map = Basemap(width=1100000,height=900000,
#         rsphere=(6378137.00,6356752.3142),
#         resolution='h',projection='lcc',
#         lat_1=45.,lat_2=55,lat_0=42,lon_0=25)


map = Basemap(projection='merc',llcrnrlat=35,urcrnrlat=45,
            llcrnrlon=23,urcrnrlon=32,lat_ts=20,resolution='h')


#TSS
# top right corner
lon0=30.0 ; lat0=42.0
# top left corner
lon1=25.0 ; lat1=42.0
# bottom left corner
lon2=25.0 ; lat2=39.0
# bottom right corner
lon3=30.0 ; lat3=39.0

#generate the new grid
lonp = np.array([lon0, lon1, lon2, lon3])
latp = np.array([lat0, lat1, lat2, lat3])

beta = np.array([1, 1, 1, 1])

print('generating the grid')
hgrd = pyroms.grid.Gridgen(lonp, latp, beta, (Mp+3,Lp+3), proj=map)
lonv, latv = list(map(hgrd.x_vert, hgrd.y_vert, inverse=True))
hgrd = pyroms.grid.CGrid_geo(lonv, latv, map)

# if you want to us ethe graphical interface
#map.drawcoastlines()
#xp, yp = map(lonp, latp)
#bry = pyroms.hgrid.BoundaryInteractor(xp, yp, beta, shp=(Mp+3,Lp+3), proj=map)
#hgrd = bry.grd


# we do not need this part untill topo
# for xx,yy in map.coastpolygons:
#     xa = np.array(xx, np.float32)
#     ya = np.array(yy,np.float32)
#     vv = np.zeros((xa.shape[0],2))
#     vv[:, 0] = xa
#     vv[:, 1] = ya
#     hgrd.mask_polygon(vv,mask_value=0)
#

# Edit the land mask interactively.
#pyroms.grid.edit_mask_mesh(hgrd, proj=map)
#edit_mask_mesh_ij is a faster version using imshow... but no map projection.
# coast = pyroms.utility.get_coast_from_map(map)
# pyroms.grid.edit_mask_mesh_ij(hgrd, coast=coast)

# datadir = '/okyanus/users/milicak/dataset/MOM6/Arctic_Copernicus_oldlonlat//'
# topo = np.loadtxt(os.path.join(datadir, 'etopo20data.gz'))
# lons = np.loadtxt(os.path.join(datadir, 'etopo20lons.gz'))
# lats = np.loadtxt(os.path.join(datadir, 'etopo20lats.gz'))

# depth positive
# topo = -topo

# fix minimum depth
# hmin = 5
# topo = np.where(topo < hmin, hmin, topo)

# interpolate new bathymetry
# lon, lat = np.meshgrid(lons, lats)
# h = griddata((lon.flat,lat.flat),topo.flat,(hgrd.lon_rho,hgrd.lat_rho), method='linear')
# h = griddata((lon.flat,lat.flat),topo.flat,(hgrd.lon_rho,hgrd.lat_rho),
#              method='nearest')

# insure that depth is always deeper than hmin
# h = np.where(h < hmin, hmin, h)

# set depth to hmin where masked
# idx = np.where(hgrd.mask_rho == 0)
# h[idx] = hmin

# save raw bathymetry
# hraw = h.copy()

# vertical coordinate
# theta_b = 2
# theta_s = 7.0
# Tcline = 50
# N = 30
# vgrd = pyroms.vgrid.s_coordinate_4(h, theta_b, theta_s, Tcline, N, hraw=hraw)
#
# ROMS grid
# grd_name = 'TSS'
# grd = pyroms.grid.ROMS_Grid(grd_name, hgrd, vgrd)

# write grid to netcdf file
# print('writing the grid')
# pyroms.grid.write_ROMS_grid(grd, filename='roms_grid_TSS.nc')

# to load the grid file back
# grd=pyroms.grid.get_ROMS_grid(grd_name,hist_file='roms_grid_Arctic_Copernicus.nc',grid_file='roms_grid_Arctic_Copernicus.nc')

# lon_rho = hgrd.lon_rho[1:-1,1:-1] # Cell centers (drop outside row and column)
# lat_rho = hgrd.lat_rho[1:-1,1:-1] # Cell centers (drop outside row and column)
# foo = xr.Dataset({'lon_rho':(['nj', 'ni'], lon_rho),'lat_rho':(['nj', 'ni'], lat_rho)})
# foo.to_netcdf('roms_grid_orig.nc')

# generate the mask
#for verts in map.coastsegs:
#    hgrd.mask_polygon(verts)
# alternate version from johan.navarro.padron

# plt.figure()
# lon, lat = map(hgrd.lon_rho,hgrd.lat_rho)
# map.fillcontinents(color='grey');
# map.plot(lon,lat,'k');
# map.plot(np.transpose(lon),np.transpose(lat),'k');
#

nj,ni = hgrd.lon_rho.shape
nj -=2; ni -=2
print('nj=%i, ni=%i'%(nj,ni))

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
lon[::2,::2] = hgrd.lon_psi[:,:] # Cell corners
lon[1::2,1::2] = hgrd.lon_rho[1:-1,1:-1] # Cell centers (drop outside row and column)
lon[1::2,::2] = hgrd.lon_u[1:-1,:] # U-points (drop outside row)
lon[::2,1::2] = hgrd.lon_v[:,1:-1] # V-points (drop outside column)
lat[::2,::2] = hgrd.lat_psi[:,:] # Cell corners
lat[1::2,1::2] = hgrd.lat_rho[1:-1,1:-1] # Cell centers (drop outside row and column)
lat[1::2,::2] = hgrd.lat_u[1:-1,:] # U-points (drop outside row)
lat[::2,1::2] = hgrd.lat_v[:,1:-1] # V-points (drop outside column)

# Approximate edge lengths as great arcs
R = 6370.e3 # Radius of sphere
dx[:,:] = R*angle_p1p2( (lat[:,1:],lon[:,1:]), (lat[:,:-1],lon[:,:-1]) )
dy[:,:] = R*angle_p1p2( (lat[1:,:],lon[1:,:]), (lat[:-1,:],lon[:-1,:]) )

# Approximate angles using centered differences in interior
lon1 = lon
cos_lat = np.cos(np.radians(lat))
angle1 = np.zeros(lat.shape)
angle2 = np.zeros(lat.shape)
# Compute it twice to recover from dateline problems, if any
angle1[:,1:-1] = np.arctan2( (lat[:,2:] - lat[:,:-2]) , ((lon[:,2:] - lon[:,:-2]) * cos_lat[:,1:-1]) )
angle1[:, 0  ] = np.arctan2( (lat[:, 1] - lat[:, 0 ]) , ((lon[:, 1] - lon[:, 0 ]) * cos_lat[:, 0  ]) )
angle1[:,-1  ] = np.arctan2( (lat[:,-1] - lat[:,-2 ]) , ((lon[:,-1] - lon[:,-2 ]) * cos_lat[:,-1  ]) )
lon = np.where(lon < 0., lon+360, lon)
angle2[:,1:-1] = np.arctan2( (lat[:,2:] - lat[:,:-2]) , ((lon[:,2:] - lon[:,:-2]) * cos_lat[:,1:-1]) )
angle2[:, 0  ] = np.arctan2( (lat[:, 1] - lat[:, 0 ]) , ((lon[:, 1] - lon[:, 0 ]) * cos_lat[:, 0  ]) )
angle2[:,-1  ] = np.arctan2( (lat[:,-1] - lat[:,-2 ]) , ((lon[:,-1] - lon[:,-2 ]) * cos_lat[:,-1  ]) )
angle = np.maximum(angle1, angle2)
lon = lon1

area = dx[:-1,:]*dy[:,:-1]
# area = area[:-1,:-1]

# Create a mosaic file
rg = scipy.io.netcdf_file('ocean_hgrid.nc','w')
# Dimensions
rg.createDimension('nx',sni)
rg.createDimension('nxp',sni+1)
rg.createDimension('ny',snj)
rg.createDimension('nyp',snj+1)
# rg.createDimension('string',255)
rg.createDimension('string',5)
# Variables
hx = rg.createVariable('x','float32',('nyp','nxp',))
hx.units = 'degrees'
hy = rg.createVariable('y','float32',('nyp','nxp',))
hy.units = 'degrees'
hdx = rg.createVariable('dx','float32',('nyp','nx',))
hdx.units = 'meters'
hdy = rg.createVariable('dy','float32',('ny','nxp',))
hdy.units = 'meters'
harea = rg.createVariable('area','float32',('ny','nx',))
harea.units = 'meters^2'
hangle = rg.createVariable('angle_dx','float32',('nyp','nxp',))
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

