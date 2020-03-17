import gcsfs                                                                    
from matplotlib import path   
                                                                                
def read_data(uri):                                                             
    gcs = gcsfs.GCSFileSystem(token='anon')                                     
    ds = xr.open_zarr(gcs.get_mapper(uri), consolidated=True)                   
    return ds                                                                   
                                                                                
                                                                                
lon_region1 = np.array([ 0.7632,    8.1791,   14.5322,   21.6440,   26.3827,     
                       24.2668, 17.4635,    8.0796,    4.0217,    0.7632]);
lat_region1 = np.array([81.0571,   82.7292,   83.5158,   83.5447,   82.5843,
                        81.2011, 80.8200,   79.2105,   79.1073,   81.0571]);



dfc = read_data('gs://cmip6/CMIP/NCAR/CESM2/historical/r1i1p1f1/Omon/thetao/gn/')
dfs = read_data('gs://cmip6/ScenarioMIP/NCAR/CESM2/ssp585/r1i1p1f1/Omon/thetao/gn/')   

lons_lats_vect = np.column_stack((lon_region1, lat_region1)) # Reshape coordinates

p = path.Path(lons_lats_vect)  

lon = np.copy(dfs.lon)  
lat = np.copy(dfs.lat)  
lons_lats_model = np.column_stack((lon.flatten(),lat.flatten())) 
mask = p.contains_points(lons_lats_model)                        
mask.shape = lon.shape
mask = np.multiply(mask, 1)     

# last 20 years of simulations 2081-2100
dfs = dfs.thetao[-240:,:,:,:]   
# years between 1981-2000
dfc = dfc.thetao[(1980-1850+1)*12:(1980-1850+1)*12+240,:,:,:]
dfs = dfs.mean('time')
dfc = dfc.mean('time')
tempc = dfc.where(mask==1)
temps = dfs.where(mask==1)

tempcmean = tempc.mean(('nlat','nlon'),skipna=True)   
tempsmean = temps.mean(('nlat','nlon'),skipna=True)   

tempcmean.to_netcdf('CESM_temp_control.nc')
tempsmean.to_netcdf('CESM_temp_ssp585.nc')

# thetao = read_data('gs://cmip6/CMIP/IPSL/IPSL-CM6A-LR/historical/r10i1p1f1/Omon/thetao/gn/')
                                                                                
                                                                                
# df = pd.read_csv('https://storage.googleapis.com/pangeo-cmip6/pangeo-cmip6-zarr-consolidated-stores.csv')
# df.head()                                                                     
                                                                                
                   #    'activity_id'+'institution_id'+'source_id'+experiment_id'+'memeber_id'+'table_id'+'vari
# able_id'+'grid_label'
# pr = read_data('gs://cmip6/ScenarioMIP/NCAR/CESM2/ssp585/r1i1p1f1/Amon/pr/gn/')
                                                                                
# df_ssp585_pr = df[(df.experiment_id == 'ssp585') & (df.variable_id == 'pr') & (df.institution_id == 'NCAR') &
#  (df.source_id == 'CESM2') & (df.table_id == 'Amon')]
                                                                                                               

