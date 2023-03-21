import numpy as np                                                              
import xarray as xr                                                             
                                                                                
df = xr.open_dataset('/work/opa/da01720/Experiments/CMEMS2_n4.2/bs-test-int-bdy4/model/domain_cfg.nc')
                                                                                
# new version                                                                   
df.bathy_metry[0,30,70:79] = 30.0                                               
# make the sill                                                                 
df.bathy_metry[0,30,74:77] = 60.0                                               
# dig after the sill                                                            
df.bathy_metry[0,31,74:77] = 67.0                                               
# extent the channel walls to east                                              
df.bathy_metry[0,31,77:79] = 42.0                                               
# make deeper region a little favor to left                                     
df.bathy_metry[0,32,77] = df.bathy_metry[0,32,78]                               
df.to_netcdf('domain_cfg_MIv3.nc')                                              
                                                                                
                                                                                
                                                                                
                                                                                
# first create a shelf                                                          
df.bathy_metry[0,30,68:79] = 30.0                                               
# make the sill                                                                 
df.bathy_metry[0,30,74:77] = 60.0                                               
# dig after the sill                                                            
df.bathy_metry[0,31,74:77] = 67.0                                               
# extent the channel walls to east                                              
df.bathy_metry[0,31,77:79] = 42.0                                               
                                                                                
df.to_netcdf('domain_cfg_MIv1.nc')                                              
                                                                                
# narrow the exit                                                               
df.bathy_metry[0,29:32,76] = df.bathy_metry[0,29:32,77]                         
df.to_netcdf('domain_cfg_MIv2.nc')                                              
                                                                                                                      
                                                                                                                      
                                                                                                                      
