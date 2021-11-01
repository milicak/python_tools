import os                                                                   
import numpy as np                                                          
import numpy.ma as ma                                                       
import glob                                                                 
import xarray as xr                                                         
import sys                                                                  
import pandas as pd                                                         
os.environ['NUMPY_EXPERIMENTAL_ARRAY_FUNCTION'] = '0'                       
                                                                            
from matplotlib import pyplot as plt                                        
import numpy as np                                                          
import xarray as xr                                                         
import gsw                                                                  
from xhistogram.xarray import histogram                                     
                                                                            
def vertical_rebin(data, bin_data, bins, dz, vert_dim="st_ocean"):          
    nanmask = np.isnan(data)                                                
    # Should we also check the bin data for nans?                           
    full_sum = histogram(                                                   
        bin_data.where(~nanmask),                                           
        bins=[bins],                                                        
        weights=(data * dz).where(~nanmask),                                
        dim=[vert_dim],                                                     
    )                                                                       
    return full_sum                                                         
                                                                            
def vertical_rebin_wrapper(                                                     
    ds,                                                                         
    bin_data_name,                                                              
    bins,                                                                       
    dz_name="dz",                                                               
    vert_dim="st_ocean",                                                        
    return_average=True,                                                        
    debug=False,                                                                
):                                                                              
    """A wrapper for the core functionality in `vertical_rebin`.                
    Accepts datasets and calculates the average over the new depth coordinates. 
    """                                                                         
    ds = ds.copy()                                                              
    ds_rebinned = xr.Dataset()                                                  
                                                                                
    ones = xr.ones_like(ds[dz_name])                                            
                                                                                
    dz_rebinned = vertical_rebin(                                               
        ones,                                                                   
        ds[bin_data_name],                                                      
        bins,                                                                   
        ds[dz_name],                                                            
        vert_dim=vert_dim,                                                      
    )                                                                           
    for var in ds.data_vars:                                                    
        ds_rebinned[var] = vertical_rebin(                                      
            ds[var], ds[bin_data_name], bins, ds[dz_name], vert_dim=vert_dim    
        )                                                                       
    if return_average:                                                          
        ds_rebinned = (                                                         
            ds_rebinned / dz_rebinned                                           
        )  # this might cause a lot of overhead...i can try to deactivate if the save fails.
                                                                                
    ds_rebinned[dz_name] = dz_rebinned                                          
                                                                                
    return ds_rebinned                                                          
                                                                                
                                                                                
gr = xr.open_dataset('/work/opa/sc02915/data/bs-simu_6.6_conf/mesh_mask_Aug19_31-5-25.nc')
gr = gr.rename({'t': 'time', 'z': 'depth', 'x': 'lon', 'y': 'lat'})

# from 1993 - 2019
root_folder = '/data/inputs/metocean/historical/model/ocean/CMCC/CMEMS/reanalysis/day/'
files = root_folder + '*/*/*RFVL*BSe3r1*'
# for 2020
root_folder = '/data/products/BSFS/bs-rea_v3.1/cmems/reanalysis_daily_mean/'
files = root_folder + np.str(2020) +'/*RFVL*BSe3r1*'      


ls1 = sorted(glob.glob(files))
root_folder = '/work/opa/mi19918/Projects/nemo/BS/MOC_data/'
# from 1993 - 2019
files = root_folder + 'sigma*'
# for 2020
files = root_folder + 'sigma2_BSe3r1_2020*'
ls2 = sorted(glob.glob(files))


                                                                                
# select bins for sigma2 or neutral                                             
bins = np.concatenate((np.array([17]),np.arange(17.5, 27, 0.1),np.array([27,28,30,32,35])))  


comp = dict(zlib=True, complevel=5)                                             
                                                                                
for ind in np.arange(0,len(ls1)):                                              
    dfs = xr.open_dataset(ls2[ind])
    dfv = xr.open_dataset(ls1[ind])
    # Mehmet
    # this was a mistake should be added the line below but I didnt
    # dfs['depth'] = -dfs.depth 
    # dfv = dfv.rename_dims({'i': 'XC', 'j': 'YC', 'i_g': 'XG', 'j_g': 'YG', 'k': 'Z', 'k_l': 'Zl', 'k_p1': 'Zp1', 'k_u': 'Zu'})
    dfv1 = dfv.isel(time=0)                                                     
    dfs['VVEL'] = dfs.sigma2                                                    
    dfs['drF'] = dfs.depth                                                          
    lon = dfs['lon'].values                                                      
    lat = dfs['lat'].values                                                      
    zlev  = dfs['depth'].values                                                     
    tmp1 = xr.DataArray(np.copy(gr.e3t_1d.isel(time=0)), coords={'depth': zlev},                    
                 dims=['depth'])                                                    
    tmp2 = xr.DataArray(np.copy(gr.e1t.isel(time=0)), coords={'lat': lat, 'lon': lon},         
             dims=['lat', 'lon'])                                                 
    me = xr.DataArray(np.copy(dfv1.vo), coords={'lat': lat, 'lon': lon,     
                                    'depth': zlev},                                 
                 dims=['depth', 'lat', 'lon'])                                        
    dfs['VVEL'] = me                                                            
    dfs['drF'] = tmp1                                                           
    df_rebinned = vertical_rebin_wrapper(dfs,                                   
                                         'sigma2',                              
                                         bins,                                  
                                         dz_name='drF',                         
                                         vert_dim='depth')                          
    df_rebinned = df_rebinned.fillna(0)                                         
    df_rebinned['dxC'] = tmp2                                                   
    voltr = df_rebinned.VVEL*df_rebinned.drF*df_rebinned.dxC                    
    voltr = voltr.to_dataset(name='vol_sigma_tr')                               
    voltr = voltr.sum('lon')                                                     
    # print(ls1[ind][79:87])
    # outname = '/work/opa/mi19918/Projects/nemo/BS/MOC_data/moc_meridional_sigma2_BSe3r1_' + ls1[ind][79:87] + '.nc'
    print(ls1[ind][65:73]) 
    outname = '/work/opa/mi19918/Projects/nemo/BS/MOC_data/moc_meridional_sigma2_BSe3r1_' + ls1[ind][65:73] + '.nc'
    encoding = {var: comp for var in voltr.data_vars}                           
    voltr.to_netcdf(outname, encoding=encoding)                                   
                                                                                
                                                                                





