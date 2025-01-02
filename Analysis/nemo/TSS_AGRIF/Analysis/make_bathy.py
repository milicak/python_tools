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

# nemo domain corners and resolution (1/100 = 0.01 degrees)
lon_nemo = np.arange(24.95,30.01,0.0075)
lat_nemo = np.arange(39.67,42.01,0.0075)
lon_nemo, lat_nemo = np.meshgrid(lon_nemo,lat_nemo)

grid_z1 = griddata(pnts, depth.flatten(), (lon_nemo, lat_nemo), method='linear')

# izmit korfezi
grid_z1[142:144,642:656]=-25.0

# open Bosphorus
grid_z1[182,540] = -60
grid_z1[182,541] = -60
grid_z1[182,542] = -60
grid_z1[183,541] = -60
grid_z1[183,542] = -60
grid_z1[184,542] = -60
grid_z1[184,543] = -60
grid_z1[185,543] = -60
grid_z1[185,544] = -60
grid_z1[186,543] = -60
grid_z1[186,544] = -60
grid_z1[187,544] = -60
grid_z1[187,545] = -60
grid_z1[187,546] = -60
grid_z1[188,546] = -60
grid_z1[189,546] = -60
grid_z1[189,547] = -60
grid_z1[190,547] = -60
grid_z1[181,541] = -60
# open Dardanelles
grid_z1[60:73,191:193] = -60.0
grid_z1[72,193] = -55.0
# prince islands
grid_z1[164,546:548] = -30.0
grid_z1[159:161,554] = -30.0

grid_z1[105,419] = -15.0

omask = np.ones(grid_z1.shape);
omask[np.isnan(grid_z1)==1]=0



# select and ocean point
new_mask = ice9(400,150,omask,False,False)
grid_z1[new_mask==0] = np.nan
# minimum and maximum depths
bb = np.where((grid_z1>-5) & (new_mask==1))
grid_z1[bb] = -5;
bb = np.where((grid_z1<-2000) & (new_mask==1))
grid_z1[bb] = -2000;


grid_z1 = -grid_z1
grid_z1[np.isnan(grid_z1)] = 0.0

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

