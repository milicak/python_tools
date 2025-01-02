import numpy as np
import xarray as xr
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

df = xr.open_dataset('ocean_hgrid_highres.nc')

lon = df.x[:,0]
lat = df.y[0,:]

# Dardanelles
lon1 = 26
lon2 = 27
lat1 = 40
lat2 = 40.5
# Bosphorus
lon3 = 28.9
lon4 = 29.2
lat3 = 40.95
lat4 = 41.25

x1 = lon.where((lon<lon1))
x1 = x1.dropna(dim='nyp')
x1 = np.copy(x1[::10])
x2 = lon.where((lon<lon2) & (lon>lon1))
x2 = x2.dropna(dim='nyp')
x2 = np.copy(x2)
x3 = lon.where((lon<lon3) & (lon>lon2))
x3 = x3.dropna(dim='nyp')
x3 = np.copy(x3[::10])
x4 = lon.where((lon<lon4) & (lon>lon3))
x4 = x4.dropna(dim='nyp')
x4 = np.copy(x4)
x5 = lon.where((lon>lon4))
x5 = x5.dropna(dim='nyp')
x5 = np.copy(x5[::10])

y1 = lat.where((lat<lat1))
y1 = y1.dropna(dim='nxp')
y1 = np.copy(y1[::10])
y2 = lat.where((lat<lat2) & (lat>lat1))
y2 = y2.dropna(dim='nxp')
y2 = np.copy(y2)
y3 = lat.where((lat<lat3) & (lat>lat2))
y3 = y3.dropna(dim='nxp')
y3 = np.copy(y3[::10])
y4 = lat.where((lat<lat4) & (lat>lat3))
y4 = y4.dropna(dim='nxp')
y4 = np.copy(y4)
y5 = lat.where((lat>lat4))
y5 = y5.dropna(dim='nxp')
y5 = np.copy(y5[::10])


lon = np.concatenate((x1,x2,x3,x4,x5))
lat = np.concatenate((y5,y4,y3,y2,y1))
dnmx,dnmy = np.meshgrid(lon,lat)
dnmx = (np.transpose(dnmx))
dnmy = (np.transpose(dnmy))
snj, sni = dnmx.shape
sni -= 1
snj -= 1
lon = np.copy(dnmx)
lat = np.copy(dnmy)
tmpx = lon
tmpy = lat

# Approximate edge lengths as great arcs
tmpdx = np.zeros((snj+1,sni))
tmpdy = np.zeros((snj,sni+1))

R = 6370.e3 # Radius of sphere
tmpdx[:,:] = R*angle_p1p2( (lat[:,1:],lon[:,1:]), (lat[:,:-1],lon[:,:-1]) )
tmpdy[:,:] = R*angle_p1p2( (lat[1:,:],lon[1:,:]), (lat[:-1,:],lon[:-1,:]) )
tmparea = tmpdx[:-1,:]*tmpdy[:,:-1]
tmpangle = np.ones((snj+1,sni+1))*-90.0

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
hx[:] = tmpx
hy[:] = tmpy
hdx[:] = tmpdx
hdy[:] = tmpdy
harea[:] = tmparea
hangle[:] = tmpangle
htile[:5] = 'tile1'
rg.close()

