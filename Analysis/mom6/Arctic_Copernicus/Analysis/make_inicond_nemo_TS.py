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
from PyCNAL_regridding import * 
                                                           

root_folder = '/data/products/OMIP/OMIP2_ORCA025/'
fname1 = 'OMIP2.025d.01_1m_19580101_19581231_grid_T.nc'
ds = xr.open_dataset(root_folder + fname1)
mom_dir = '/work/opa/mi19918/Projects/mom6/Arctic_Copernicus/INPUT/'
path_regional_grid = mom_dir + './ocean_hgrid.nc'                     

variable = 'thetao'

ds[variable] = ds.thetao.where(ds.thetao!=0)
df = flood_kara(ds[variable][0,:,:,:], xdim='x_grid_T', ydim='y_grid_T', zdim='deptht') 
df = df.to_dataset(name='temp')  
df['nav_lon_grid_T'] = ds['nav_lon_grid_T']  
df['nav_lat_grid_T'] = ds['nav_lat_grid_T']  
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

# for regridding                                                        
df = df.rename({'deptht':'depth'})
temp = np.copy(df.temp)                                                 
dstemp = xr.DataArray(temp, coords=[df.time, df.depth, df.y, df.x],     
                            dims=["time","depth","y","x"])              
dstemp = dstemp.ffill(dim='depth') 
df['temp'] = dstemp                                                     
                                                                        
# build regridder                                                       
regridder = xe.Regridder(df, ds2, 'nearest_s2d')
                                                                        
#apply regridder                                                        
dr_out = regridder(df['temp'])                                          
# create 3d ssh variable                   
# not sure if i need to multiply with mask 
var = np.copy(dr_out)                      
nk = var.shape[1]                          


# salinity                                                              
variable = 'so'                                                    
ds[variable] = ds.so.where(ds.so!=0)
df = flood_kara(ds[variable][0,:,:,:], xdim='x_grid_T', ydim='y_grid_T', zdim='deptht') 
df = df.to_dataset(name='salt')                                         
df['nav_lon_grid_T'] = ds['nav_lon_grid_T']  
df['nav_lat_grid_T'] = ds['nav_lat_grid_T']  
df = df.rename({'nav_lon_grid_T': 'lon', 'nav_lat_grid_T': 'lat', 'y_grid_T':
                'y', 'x_grid_T': 'x'})            
df = df.rename({'deptht':'depth'})
# for regridding                                                        
salt = np.copy(df.salt)                                                 
dssalt = xr.DataArray(salt, coords=[df.time, df.depth, df.y, df.x],     
                            dims=["time","depth","y","x"])              
dssalt = dssalt.ffill(dim='depth') 
df['salt'] = dssalt                                                     
#apply regridder                                                        
dr_out = regridder(df['salt'])                                          
varsalt = np.copy(dr_out)                                               
                                                                        
time = 17.5                                                                
                                                                           
# Create a mosaic file      
fout  = mom_dir + 'NEMO_TS_IC.nc'
rg = scipy.io.netcdf_file(fout,'w')                            
# Dimensions                                                               
rg.createDimension('time', None)                                           
rg.createDimension('z_l',nk)                                               
rg.createDimension('nxp',ni)                                               
rg.createDimension('nyp',nj)                                               
# Variables                                                                
hnx = rg.createVariable('nxp', 'int32', ('nxp',))                          
hny = rg.createVariable('nyp', 'int32', ('nyp',))                          
hx = rg.createVariable('lon','float32',('nyp','nxp',))                     
hx.units = 'degrees'                                                       
hy = rg.createVariable('lat','float32',('nyp','nxp',))                     
hy.units = 'degrees'                                                       
hz = rg.createVariable('z_l','float32',('z_l',))                           
hz.units = 'meters'                                                        
tempvar  = rg.createVariable('temp','float32',('time','z_l','nyp','nxp',)) 
tempvar.units = 'celcius'                                                  
tempvar.missing_val = 1e20                                                 
tempvar._FillValue = 1e20                                                  
saltvar  = rg.createVariable('salt','float32',('time','z_l','nyp','nxp',)) 
saltvar.units = 'psu'                                                      
saltvar.missing_val = 1e20                                                 
saltvar._FillValue = 1e20                                                  
htime = rg.createVariable('time','float32',('time',))                      
# htime = rg.createVariable('time', 'int', ('time'))                       
htime.units = 'days since 1958-01-01 00:00:00'                             
# Values                                                                   
hx[:] = lon_rho                                                            
hy[:] = lat_rho                                                            
hz[:] = np.copy(df.depth)                                                  
tempvar[:] = var                                                           
saltvar[:] = varsalt                                                       
hnx[:] = np.arange(0,ni)                                                   
hny[:] = np.arange(0,nj)                                                   
htime = time                                                               
rg.close()                                                                 
                                                                           
