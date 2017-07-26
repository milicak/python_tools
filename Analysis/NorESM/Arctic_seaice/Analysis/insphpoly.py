import numpy as np
from matplotlib.path import Path


# plot the map and get the points
# middle button to turn it off
vertices=plt.ginput()
vertices=np.copy(vertices)
vertices=np.append(vertices,vertices[0,:])
vertices=np.reshape(vertices,(vertices.size/2,2)
                    # now I have lon lat values for the polygon
lon = nc_read(grdname,'plon')
lat = nc_read(grdname,'plat')

# or you can create points
#lon, lat = np.meshgrid(x, y)  # X, Y are 2D ndarrays

lonlat = np.dstack((lon, lat))
lonlat_flat = lonlat.reshape((-1, 2))

mpath = Path( vertices ) # the vertices of the polygon
mask_flat = mpath.contains_points(lonlat_flat)
mask = mask_flat.reshape(lon.shape)
