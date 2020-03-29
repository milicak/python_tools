import xarray as xr                                                      
import numpy as np
                                                                         
root_folder = '/shared/projects/uniklima/globclim/milicak/dataset/ERA5/'
                                                                         
# varname = 'era5_global_2m_temperature_'                                
# var = 't2m'                                                            
# var2 = 't2m'                                                           
# coeff = 1.0                                                            
# varname = 'era5_global_10m_u_component_of_wind_'                       
# var = 'u10'                                                            
# var2 = 'u10'                                                           
# coeff = 1.0                                                            
# varname = 'era5_global_10m_v_component_of_wind_'                       
# var = 'v10'                                                            
# var2 = 'v10'                                                           
# coeff = 1.0                                                            
# varname = 'era5_specific_humidity_'                                    
# var = 'q'                                                              
# var2 = 'spfh2m'                                                        
# coeff = 1.0                                                            
varname = 'era5_global_surface_thermal_radiation_downwards_'           
var = 'strd'                                                           
var2 = 'dlw'                                                           
coeff = 3600.0                                                         
# varname = 'era5_global_surface_solar_radiation_downwards_'               
# var = 'ssrd'                                                             
# var2 = 'dsw'                                                             
# coeff = 3600.0                                                           


for year in range(2003,2004):                                  
    fname = root_folder + varname + str(year) + '.nc'          
    df = xr.open_dataset(fname)                                
    ds = df.resample(time='6H')
    dsm = ds.mean('time')
    dnm = np.flip(np.copy(dsm[var])/coeff,axis=1)               
    outputfile = root_folder + 'ERA5_' + var2 + '_' + str(year)
    print(outputfile)                                          
    dnm.astype('>f4').tofile(outputfile)                       
                                                               
