import numpy as np                                                              
import scipy.io                                                                 
                                                                                
                                                                                
df = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/MOM6/NWA25/ocean_topog_v2.nc')
                                                                                

offset = 1152
hraw = np.copy(df.depth[:offset,:])
nj, ni = hraw.shape                                                             
                                                                                
# Create a topography file                                                      
rg = scipy.io.netcdf_file('ocean_topog.nc','w')                                 
# Dimensions                                                                    
rg.createDimension('nx',ni)                                                     
rg.createDimension('ny',nj)                                                     
# Variables                                                                     
hdepth = rg.createVariable('depth','float32',('ny','nx',))                      
hdepth.units = 'm'                                                              
# Values                                                                        
hdepth[:] = hraw #[0,1:-1,1:-1]                                                 
rg.close()                                                                      
                                                                                
cmnd1 ="ncap2 -s " + " 'defdim(" + '"ntiles"' + ",1)' " + " ocean_topog.nc dnm.nc"
os.system(cmnd1)                                                                
cmnd1 = 'mv dnm.nc ocean_topog.nc'                                              
os.system(cmnd1)                                                                
cmnd1 = 'cp ocean_topog.nc /Volumes/A1/workdir/milicak/datasets/MOM6/NWA25_v2/.'                  
os.system(cmnd1)                                                                
                                                                                                                 
