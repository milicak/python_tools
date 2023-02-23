import os                                                                 
import numpy as np                                                        
import xesmf as xe
import xarray as xr                                                       
import scipy.io                                                           
import netCDF4                                                            
import matplotlib.pyplot as plt                                           
                                                                          
                                                                          
root_folder = '/Volumes/A1/workdir/milicak/datasets/MOM6/NWA25/'
mom_dir = '/Volumes/A1/workdir/milicak/datasets/MOM6/NWA25/'
path_regional_grid = mom_dir + './ocean_hgrid.nc'                         
                                                                          
df2 = xr.open_dataset(path_regional_grid)                                 
lon_rho = np.copy(df2['x'][1::2,1::2])                                    
lat_rho = np.copy(df2['y'][1::2,1::2])                                    
nj,ni = lon_rho.shape                                                     
ds2 = xr.Dataset()
ds2['lon'] = df2['x'][1::2,1::2]
ds2['lat'] = df2['y'][1::2,1::2]

df = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/glorys/GLORYS_REANALYSIS_1996-01-01.nc')
df2['lon'] = df['longitude']
df2['lat'] = df['latitude']

regridder = xe.Regridder(df2, ds2, 'bilinear', reuse_weights=True)
#apply regridder                   
dr_out = regridder(df[zos])   
dr_out = dr_out.fillna(0)          
                                                                          
ssh = np.copy(dr_out)                                                     
                                                                          
time = 17.5                                                               

# Create a mosaic file                                                          
fout  = mom_dir + 'glorys_ssh_1996_IC.nc'                                             
rg = scipy.io.netcdf_file(fout,'w')                                             
# Dimensions                                                                    
rg.createDimension('time', None)                                                
rg.createDimension('nxp',ni)                                                    
rg.createDimension('nyp',nj)                                                    
# Variables                                                                     
hnx = rg.createVariable('nxp', 'int32', ('nxp',))                               
hny = rg.createVariable('nyp', 'int32', ('nyp',))                               
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
htime.units = 'days since 1996-01-01 00:00:00'                                  
# Values                                                                        
hx[:] = lon_rho                                                                 
hy[:] = lat_rho                                                                 
sshvar[:] = ssh                                                                 
hnx[:] = np.arange(0,ni)                                                        
hny[:] = np.arange(0,nj)                                                        
htime = time                                                                    
rg.close()                                                                      

