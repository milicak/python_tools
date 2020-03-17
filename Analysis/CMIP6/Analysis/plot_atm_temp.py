import gcsfs                                                                    
                                                                                
def read_data(uri):                                                             
    gcs = gcsfs.GCSFileSystem(token='anon')                                     
    ds = xr.open_zarr(gcs.get_mapper(uri), consolidated=True)                   
    return ds                                                                   
                                                                                
                                                                                
                   #    'activity_id'+'institution_id'+'source_id'+experiment_id'+'memeber_id'+'table_id'+'variable_id'+'grid_label'
df = read_data('gs://cmip6/ScenarioMIP/NCAR/CESM2/ssp585/r1i1p1f1/Amon/tas/gn/') 
ds_time = pd.date_range('2015-01-01', freq='M', periods=12 * 86)
df['time'] = ds_time
df1 = df.sel(time=slice('2081-01-01','2100-12-31')) # 20 years - rcp

dc = read_data('gs://cmip6/CMIP/NCAR/CESM2/historical/r1i1p1f1/Amon/tas/gn/') 
dh_time = pd.date_range('1850-01-01', freq='M', periods=12 * 165)
dc['time'] = dh_time
dc1 = dc.sel(time=slice('1981-01-01','2000-12-31'))

dc1 = dc1.mean('time')
df1 = df1.mean('time')
