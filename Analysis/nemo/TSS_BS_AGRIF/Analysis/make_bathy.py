import numpy as np
import scipy.io
from scipy.interpolate import griddata
import xarray as xr


def ice9(i, j, source, xcyclic=True, tripolar=True):
    """
    An iterative (stack based) implementation of "Ice 9".
    The flood fill starts at [j,i] and treats any positive value of "source" as
    passable. Zero and negative values block flooding.
    xcyclic = True allows cyclic behavior in the last index. (default)
    tripolar = True allows a fold across the top-most edge. (default)
    Returns an array of 0's and 1's.
    """
    wetMask = 0 * source
    (nj, ni) = wetMask.shape
    stack = set()
    stack.add((j, i))
    while stack:
        (j, i) = stack.pop()
        if wetMask[j, i] or source[j, i] <= 0:
            continue
        wetMask[j, i] = 1
        if i > 0:
            stack.add((j, i - 1))
        elif xcyclic:
            stack.add((j, ni - 1))
        if i < ni - 1:
            stack.add((j, i + 1))
        elif xcyclic:
            stack.add((j, 0))
        if j > 0:
            stack.add((j - 1, i))
        if j < nj - 1:
            stack.add((j + 1, i))
        elif tripolar:
            stack.add((j, ni - 1 - i))  # Tri-polar fold
    return wetMask


ibcaofile = '/okyanus/users/milicak/dataset/world_bathy/GEBCO_2014_1D.nc';

df = xr.open_dataset(ibcaofile)
depth = np.reshape(np.copy(df.z),(21600,43200))
depth = np.flipud(depth)

x_range = np.copy(df.x_range)
y_range = np.copy(df.y_range)
lon = np.linspace(x_range[0],x_range[1],43200)
lat = np.linspace(y_range[0],y_range[1],21600)

# cut our region
lat = lat[15450:16600]
lon = lon[24300:26700]
depth = depth[15450:16600,24300:26700]
depth[depth>0] = np.nan
lon_gebco, lat_gebco = np.meshgrid(lon,lat)
pnts = np.transpose(np.array((lon_gebco.flatten(),lat_gebco.flatten())))

# nemo domain corners and resolution (1/50 = 0.02 degrees)
lon_nemo = np.arange(22.5,42.01,0.02)
lat_nemo = np.arange(38.8,47.51,0.02)
lon_nemo, lat_nemo = np.meshgrid(lon_nemo,lat_nemo)

grid_z1 = griddata(pnts, depth.flatten(), (lon_nemo, lat_nemo), method='linear')

# open Bosphorus
grid_z1[111,326] = -60
grid_z1[112,326] = -60
grid_z1[112,327] = -60
grid_z1[113,327] = -60
grid_z1[114,328] = -60
grid_z1[115,328] = -60
grid_z1[113,328] = -60
grid_z1[118,329] = -60
grid_z1[119,330] = -60
# open Dardanelles
grid_z1[73,203] = -60
grid_z1[74,203] = -60
grid_z1[75,203] = -60
grid_z1[70,197] = -60
grid_z1[70,196] = -60
grid_z1[69,196] = -60
grid_z1[68,196] = -60
grid_z1[69,195] = -60
grid_z1[68,195] = -60
grid_z1[67,195] = -60
grid_z1[67,194] = -60
grid_z1[65,193] = -60
# prince islands
grid_z1[105,332] = -30
grid_z1[104,331] = -30
grid_z1[104,330] = -30
grid_z1[105,327] = -30

omask = np.ones(grid_z1.shape);
omask[np.isnan(grid_z1)==1]=0



# select and ocean point
new_mask = ice9(400,200,omask,False,False)
grid_z1[new_mask==0] = np.nan
# minimum and maximum depths
bb = np.where((grid_z1>-5) & (new_mask==1))
grid_z1[bb] = -5;
bb = np.where((grid_z1<-2500) & (new_mask==1))
grid_z1[bb] = -2500;

ny,nx = grid_z1.shape

# Create a mosaic file
rg = scipy.io.netcdf_file('bathy_meter.nc','w')
# Dimensions
rg.createDimension('x',nx)
rg.createDimension('y',ny)
# Variables
hx = rg.createVariable('nav_lon','float32',('y','x',))
hx.units = 'degrees_east'
hx.longname = 'Longitude'
hy = rg.createVariable('nav_lat','float32',('y','x',))
hy.units = 'degrees_north'
hy.longname = 'Latitude'
hdx = rg.createVariable('Bathymetry','float32',('y','x',))
hdx.units = 'm'
hdx.longname = 'bathymetry'
# Values
hx[:] = lon_nemo
hy[:] = lat_nemo
hdx[:] = grid_z1
rg.close()
