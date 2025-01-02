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

df = xr.open_dataset('ocean_hgrid_coarse.nc')

ny, nx = df.area.shape
sni = nx*4
snj = ny*4

# lon = np.arange(0,nx-1+0.25,0.25)
# lat = np.arange(0,ny-1+0.25,0.25)

lon1 = np.zeros(sni+1)
lat1 = np.zeros(snj+1)
lon1[1:] = np.linspace(0.25,nx,sni,endpoint=True)
lat1[1:] = np.linspace(0.25,ny,snj,endpoint=True)

tmpx = df.x.interp(nxp=lon1,nyp=lat1)
tmpy = df.y.interp(nxp=lon1,nyp=lat1)

# tmpdx = df.dx.interp(nx=lon,nyp=lat1)/4
# tmpdy = df.dy.interp(nxp=lon1,ny=lat)/4

tmpangle = np.ones((snj+1,sni+1))*-90.0

# Approximate edge lengths as great arcs
tmpdx = np.zeros((snj+1,sni))
tmpdy = np.zeros((snj,sni+1))
lon = np.copy(tmpx)
lat = np.copy(tmpy)

R = 6370.e3 # Radius of sphere
tmpdx[:,:] = R*angle_p1p2( (lat[:,1:],lon[:,1:]), (lat[:,:-1],lon[:,:-1]) )
tmpdy[:,:] = R*angle_p1p2( (lat[1:,:],lon[1:,:]), (lat[:-1,:],lon[:-1,:]) )
tmparea = tmpdx[:-1,:]*tmpdy[:,:-1]

# Create a mosaic file
rg = scipy.io.netcdf_file('ocean_hgrid_highres.nc','w',version=2)
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

