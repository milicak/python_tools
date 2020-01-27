import gcsfs

def read_data(uri):
    gcs = gcsfs.GCSFileSystem(token='anon')
    ds = xr.open_zarr(gcs.get_mapper(uri), consolidated=True)
    return ds


thetao = read_data('gs://cmip6/CMIP/IPSL/IPSL-CM6A-LR/historical/r10i1p1f1/Omon/thetao/gn/')


df = pd.read_csv('https://storage.googleapis.com/pangeo-cmip6/pangeo-cmip6-zarr-consolidated-stores.csv')
df.head()

                   #    'activity_id'+'institution_id'+'source_id'+experiment_id'+'memeber_id'+'table_id'+'variable_id'+'grid_label'
pr = read_data('gs://cmip6/ScenarioMIP/NCAR/CESM2/ssp585/r1i1p1f1/Amon/pr/gn/')

df_ssp585_pr = df[(df.experiment_id == 'ssp585') & (df.variable_id == 'pr') & (df.institution_id == 'NCAR') & (df.source_id == 'CESM2') & (df.table_id == 'Amon')]
