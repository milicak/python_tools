import numpy as np
import scipy.io
import xesmf as xe


# df = xr.open_dataset('~/dataset/MOM6/NWA25_v2/ocean_topog.nc')
df = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/NWA25_v2/ocean_mask.nc')
df2 = xr.open_dataset('/okyanus/users/milicak/dataset/MOM6/NWA25_v2/ocean_hgrid.nc')
# lon_rho = np.copy(df2['x'][1::2,1::2])
# lat_rho = np.copy(df2['y'][1::2,1::2])
ds2 = df2['x'][1::2,1::2]
ds2 = ds2.to_dataset(name='lon')
ds2['lat']=df2['y'][1::2,1::2]
ds2 = ds2.rename_dims({'nxp': 'x','nyp': 'y'})

# Read and interpolate TOPO file into the grid
gr = xr.open_dataset('/okyanus/users/milicak/dataset/world_bathy/GEBCO_2014_1D.nc')
lon_gebco = np.linspace(np.copy(gr.x_range[0]),np.copy(gr.x_range[1]),np.copy(gr.dimension[0]),endpoint=True);
lat_gebco = np.linspace(np.copy(gr.y_range[0]),np.copy(gr.y_range[1]),np.copy(gr.dimension[1]),endpoint=True);
depth_gebco = np.reshape(np.copy(gr.z), (np.copy(gr.dimension[1]), np.copy(gr.dimension[0])))
lat_gebco = lat_gebco[11000:18000]
lon_gebco = lon_gebco[9000:18000]
depth_gebco = np.flipud(depth_gebco)
depth_gebco = depth_gebco[11000:18000,9000:18000]

# Create xarray dataset for depth_gebco
foo = xr.Dataset({'depth_g':(['lat','lon'],  depth_gebco)})
foo = foo.assign_coords(lat=lat_gebco,lon=lon_gebco)

# dd = foo.interp(lat=df.lat_rho,lon=df.lon_rho,method='linear')

# build regridder
regridder = xe.Regridder(foo, ds2, 'bilinear', reuse_weights=True)

#apply regridder
dr_out = regridder(foo.depth_g)
dfs = dr_out.to_dataset(name='depth')

aa = np.copy(dfs.depth.where(dfs.depth<0,0))
aa[(aa>-5) & (aa<0)]=-5
aa[aa<-6500]=-6500

aa = -aa*np.copy(df.mask)

hraw = np.copy(aa)

nj, ni = hraw.shape

# Create a topography file
rg = scipy.io.netcdf_file('/okyanus/users/milicak/dataset/MOM6/NWA25_v2/ocean_topog_gebco.nc','w')
# Dimensions
rg.createDimension('nx',ni)
rg.createDimension('ny',nj)
rg.createDimension('ntiles',1)
# Variables
hdepth = rg.createVariable('depth','float32',('ny','nx',))
hdepth.units = 'm'
# htile = rg.createVariable('tile','c',('string',))
# htile[:5] = 'tile1'
# Values
hdepth[:] = hraw #[0,1:-1,1:-1]
rg.close()

