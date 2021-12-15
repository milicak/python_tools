import numpy as np


root_folder = '/data/products/OMIP/OMIP2_ORCA025/'
fname1 = 'OMIP2.025d.01_1m_19580101_19581231_grid_T.nc'
df = xr.open_dataset(root_folder + fname1)
mom_dir = '/work/opa/mi19918/Projects/mom6/Arctic_Copernicus/INPUT/'
path_regional_grid = mom_dir + './ocean_hgrid.nc'


ds = xr.Dataset()                                                                                                                                        
ds['temp']=df["thetao"][:,:,:-2,1:-1]  
ds['salt']=df["so"][:,:,:-2,1:-1]       
ds['ssh']=df["zos"][:,:-2,1:-1]

ds = ds.rename({'nav_lon_grid_T':'lon','nav_lat_grid_T':'lat','x_grid_T':'x','y_grid_T':'y'})  
ds = ds.rename({'time_counter': 'time'})

ds.to_netcdf(mom_dir+fname1)
