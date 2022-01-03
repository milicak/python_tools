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
                                                                                
variable = 'aice'                                                               
variable1 = 'hi'                                                                
                                                                                

ds = xr.open_dataset('/data/products/OMIP/OMIP2_ORCA025/OMIP2.025d.01.cice.h.1958-01.nc')
mom_dir = '/work/opa/mi19918/Projects/mom6/Arctic_Copernicus/INPUT/'


# ice concentration should be between 0 and 1                                   
ds[variable] = ds[variable]/100.0
df = ds[variable]                                                              
df = df.to_dataset(name=variable)                                               
df['hi'] = ds[variable1]                                                       
df = df.rename({'TLON': 'lon', 'TLAT': 'lat', 'nj': 'y', 'ni': 'x'})            
                                                                                
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
regridder = xe.Regridder(df, ds2, 'nearest_s2d', reuse_weights=True)            
                                                                                
#apply regridder                          
dr_out = regridder(df[variable])          
dr_out = dr_out.fillna(0)                 
dr_out2 = regridder(df[variable1])        
dr_out2 = dr_out2.fillna(0)               
aice = np.copy(dr_out)                    
hice = np.copy(dr_out2)                   
                                          
time = 17.5                               
                                          
# Create a mosaic file                                                    
fout = mom_dir + 'NEMO_seaice_IC.nc'
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
aicevar  = rg.createVariable('aice','float32',('time','nyp','nxp',))      
aicevar.units = 'concentration 0-1'                                       
aicevar.missing_val = 1e20                                                
aicevar._FillValue = 1e20                                                 
hicevar  = rg.createVariable('hice','float32',('time','nyp','nxp',))      
hicevar.units = 'meters'                                                  
hicevar.missing_val = 1e20                                                
hicevar._FillValue = 1e20                                                 
htime = rg.createVariable('time','float32',('time',))                     
# htime = rg.createVariable('time', 'int', ('time'))                      
htime.units = 'days since 1958-01-01 00:00:00'                            
# Values                                                                  
hx[:] = lon_rho    
hy[:] = lat_rho    
aicevar[:] = aice  
hicevar[:] = hice  
htime = time       
rg.close()         

