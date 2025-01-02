import numpy as np
import scipy.io


df = xr.open_dataset('~/dataset/MOM6/Marguerite_Bay/margi_v04.nc')

hraw = np.copy(df.h[1:-1,1:-1]*df.mask_rho[1:-1,1:-1])

nj, ni = hraw.shape

# Create a topography file
rg = scipy.io.netcdf_file('/okyanus/users/milicak/dataset/MOM6/Marguerite_Bay/ocean_topog.nc','w')
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

