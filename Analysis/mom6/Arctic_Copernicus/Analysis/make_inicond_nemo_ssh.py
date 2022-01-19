import os                                                                       
import numpy as np                                                              
import xesmf as xe                                                              
import xarray as xr                                                             
import scipy.io                                                                 
from scipy.io import savemat                                                    
from scipy.io import loadmat                                                    
import matplotlib.colors as colors                                              
from scipy.signal import medfilt2d                                              
import netCDF4                                                                  
import matplotlib.pyplot as plt                                                 
from scipy.interpolate import griddata                                          
from matplotlib.path import Path                                                
#for interpolation                                                              
from scipy.spatial import cKDTree                                               
from HCtFlood.kara import flood_kara                                            
                                                                                
variable = 'zos'                                                               

root_folder = '/data/products/OMIP/OMIP2_ORCA025/'
fname1 = 'OMIP2.025d.01_1m_19580101_19581231_grid_T.nc'
ds = xr.open_dataset(root_folder + fname1)
mom_dir = '/work/opa/mi19918/Projects/mom6/Arctic_Copernicus/INPUT/'


df = ds[variable][0,:,:]
df = df.to_dataset(name=variable)                                               
df = df.rename({'nav_lon_grid_T': 'lon', 'nav_lat_grid_T': 'lat', 'y_grid_T':
                'y', 'x_grid_T': 'x'})            
                                                                                
df2 = xr.open_dataset('/work/opa/mi19918/Projects/mom6/Arctic_Copernicus/INPUT/ocean_hgrid.nc')
lon_rho = np.copy(df2['x'][1::2,1::2])                                          
lat_rho = np.copy(df2['y'][1::2,1::2])                                          
nj,ni = lon_rho.shape                                                           
ds2 = df2['x'][1::2,1::2]                                                       
ds2 = ds2.to_dataset(name='lon')                                                
ds2['lat']=df2['y'][1::2,1::2]                                                  
ds2 = ds2.rename_dims({'nxp': 'x','nyp': 'y'})                                  
                                                                                
# build regridder                                                               
regridder = xe.Regridder(df, ds2, 'nearest_s2d')                                
                                                                                
#apply regridder                          
dr_out = regridder(df[variable])          
ssh = np.copy(dr_out)
ssh =np.reshape(ssh,(1,nj,ni))
                                          
time = 17.5                               
                                          
# Create a mosaic file                                                    
fout = mom_dir + 'NEMO_ssh_IC.nc'
rg = scipy.io.netcdf_file(fout,'w')                             
# Dimensions                                                              
rg.createDimension('time', None)                                     
rg.createDimension('nxp',ni)                                         
rg.createDimension('nyp',nj)                                         
# Variables                                                          
hx = rg.createVariable('lon','float32',('nyp','nxp',))               
hx.units = 'degrees'                                                 
hy = rg.createVariable('lat','float32',('nyp','nxp',))               
hy.units = 'degrees'                                                 
sshvar  = rg.createVariable('ssh','float32',('time','nyp','nxp',))   
sshvar.units = 'meters'                                              
sshvar.missing_val = 1e20                                            
sshvar._FillValue = 1e20                                             
htime = rg.createVariable('time','float32',('time',))                
# htime = rg.createVariable('time', 'int', ('time'))                 
htime.units = 'days since 1958-01-01 00:00:00'                       
# Values                                                             
hx[:] = lon_rho                                                      
hy[:] = lat_rho                                                      
sshvar[:] = ssh                                                      
htime = time                                                         
rg.close()                                                           

