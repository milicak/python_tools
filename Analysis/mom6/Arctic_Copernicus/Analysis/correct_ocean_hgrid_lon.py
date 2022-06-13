import numpy as np
import scipy.io

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

def angle_between(v1, v2, v3):
    """Returns angle v2-v1-v3 i.e betweeen v1-v2 and v1-v3."""
    # vector product between v1 and v2
    px = v1[1] * v2[2] - v1[2] * v2[1]
    py = v1[2] * v2[0] - v1[0] * v2[2]
    pz = v1[0] * v2[1] - v1[1] * v2[0]
    # vector product between v1 and v3
    qx = v1[1] * v3[2] - v1[2] * v3[1]
    qy = v1[2] * v3[0] - v1[0] * v3[2]
    qz = v1[0] * v3[1] - v1[1] * v3[0]

    ddd = (px * px + py * py + pz * pz) * (qx * qx + qy * qy + qz * qz)
    ddd = (px * qx + py * qy + pz * qz) / np.sqrt(ddd)
    angle = np.arccos( ddd )
    return angle

def quad_area(lat, lon):
    """Returns area of spherical quad (bounded by great arcs)."""
    # x,y,z are 3D coordinates
    d2r = np.deg2rad(1.)
    x = np.cos(d2r * lat) * np.cos(d2r * lon)
    y = np.cos(d2r * lat) * np.sin(d2r * lon)
    z = np.sin(d2r * lat)
    c0 = (x[ :-1, :-1], y[ :-1, :-1], z[ :-1, :-1])
    c1 = (x[ :-1,1:  ], y[ :-1,1:  ], z[ :-1,1:  ])
    c2 = (x[1:  ,1:  ], y[1:  ,1:  ], z[1:  ,1:  ])
    c3 = (x[1:  , :-1], y[1:  , :-1], z[1:  , :-1])
    a0 = angle_between(c1, c0, c2)
    a1 = angle_between(c2, c1, c3)
    a2 = angle_between(c3, c2, c0)
    a3 = angle_between(c0, c3, c1)
    return a0 + a1 + a2 + a3 - 2. * np.pi


df = xr.open_dataset('~/dataset/MOM6/Arctic_Copernicus/ocean_hgrid.nc')

lon = np.copy(df.x)
lon1 = lon
lon[lon<0]=lon[lon<0]+360
lat = np.copy(df.y)

lon = np.transpose(lon)
lon1 = np.transpose(lon1)
lat = np.transpose(lat)

snj, sni = lon.shape

snj -= 1
sni -= 1

area = np.zeros((snj,sni))
dx = np.zeros((snj+1,sni))
dy = np.zeros((snj,sni+1))
angle = np.zeros((snj+1,sni+1))

# Approximate edge lengths as great arcs
R = 6370.e3 # Radius of sphere
dx[:,:] = R*angle_p1p2( (lat[:,1:],lon[:,1:]), (lat[:,:-1],lon[:,:-1]) )
dy[:,:] = R*angle_p1p2( (lat[1:,:],lon[1:,:]), (lat[:-1,:],lon[:-1,:]) )

# Approximate angles using centered differences in interior
cos_lat = np.cos(np.radians(lat))
angle1 = np.zeros(lat.shape)
angle2 = np.zeros(lat.shape)
# Compute it twice to recover from dateline problems, if any
angle1[:,1:-1] = np.arctan2( (lat[:,2:] - lat[:,:-2]) , ((lon[:,2:] - lon[:,:-2]) * cos_lat[:,1:-1]) )
angle1[:, 0  ] = np.arctan2( (lat[:, 1] - lat[:, 0 ]) , ((lon[:, 1] - lon[:, 0 ]) * cos_lat[:, 0  ]) )
angle1[:,-1  ] = np.arctan2( (lat[:,-1] - lat[:,-2 ]) , ((lon[:,-1] - lon[:,-2 ]) * cos_lat[:,-1  ]) )
lon1 = np.where(lon1 < 0., lon1+360, lon1)
angle2[:,1:-1] = np.arctan2( (lat[:,2:] - lat[:,:-2]) , ((lon1[:,2:] - lon1[:,:-2]) * cos_lat[:,1:-1]) )
angle2[:, 0  ] = np.arctan2( (lat[:, 1] - lat[:, 0 ]) , ((lon1[:, 1] - lon1[:, 0 ]) * cos_lat[:, 0  ]) )
angle2[:,-1  ] = np.arctan2( (lat[:,-1] - lat[:,-2 ]) , ((lon1[:,-1] - lon1[:,-2 ]) * cos_lat[:,-1  ]) )
angle = np.maximum(angle1, angle2)

area = dx[:-1,:]*dy[:,:-1]
# area = R * R * quad_area(lat, lon)

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

