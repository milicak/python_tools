import numpy as np
import scipy.io                                            

mom_dir = '/work/opa/mi19918/Projects/mom6/Arctic_Copernicus/INPUT/'
df = xr.open_dataset('/work/opa/mi19918/Projects/mom6/Arctic_Copernicus/INPUT/NEMO_TS_IC.nc')
time = 17.5                                                                

var = np.copy(df.temp)
varsalt = np.copy(df.salt)

# uniform salt
varsalt[:,:,:,:] = 35.0
templin = np.linspace(1,-1,50)
templ = np.tile(templin,(1,1750,2000,1))
var = templ.swapaxes(1,3).swapaxes(2,3)

# Create a mosaic file      
nk = df.z_l.shape[0]
ni = df.nxp.shape[0]
nj = df.nyp.shape[0]
lon_rho = np.copy(df.lon)
lat_rho = np.copy(df.lat)

fout  = mom_dir + 'NEMO_TS_IC_v3.nc'
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
hz[:] = np.copy(df.z_l)                                                  
tempvar[:] = var                                                           
saltvar[:] = varsalt                                                       
hnx[:] = np.arange(0,ni)                                                   
hny[:] = np.arange(0,nj)                                                   
htime = time                                                               
rg.close()                                                                 
