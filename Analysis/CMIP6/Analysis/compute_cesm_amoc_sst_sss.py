import gcsfs
from matplotlib import path

def read_data(uri):
    gcs = gcsfs.GCSFileSystem(token='anon')
    ds = xr.open_zarr(gcs.get_mapper(uri), consolidated=True)
    return ds


def areamean(var,area, lonvar='nlon', latvar='nlat'):
    ds = var.sum((lonvar,latvar),skipna=True)/area.sum((lonvar,latvar),skipna=True)
    return ds

def crmask(lonlat_polygon, lonlat_model, lonmodel):
    p = path.Path(lonlat_polygon)
    mask = p.contains_points(lonlat_model)
    mask.shape = lonmodel.shape
    mask = np.multiply(mask, 1)
    return mask


def areawgtmean(df, gr, mask):
    area = gr.areacello.where(mask==1)
    sst = df.tos.where(mask==1)
    sstmean = sst.groupby('time.year').mean('time')
    sstmean = sstmean*gr.areacello
    return sstmean

def areawgtmeansalt(df, gr, mask):
    area = gr.areacello.where(mask==1)
    sss = df.so[:,0,:,:].where(mask==1)
    sssmean = sss.groupby('time.year').mean('time')
    sssmean = sssmean*area
    # sssmean = sssmean*gr.areacello
    return sssmean



# Subpolar gyre area
sparea = np.genfromtxt('subpolar_gyre_lonlat.txt',usecols=(0,1))
# Subtropical gyre area
starea = np.genfromtxt('subtropical_gyre_lonlat.txt',usecols=(0,1))

dfs = read_data('gs://cmip6/ScenarioMIP/NCAR/CESM2/ssp585/r1i1p1f1/Omon/tos/gn/')
gr = xr.open_dataset('areacello_Ofx_CESM2_ssp585_r2i1p1f1_gn.nc')

landmask=np.copy(dfs.tos[0,:,:])
landmask = np.multiply(np.isnan(landmask), 1)

lon = np.copy(dfs.lon)
lat = np.copy(dfs.lat)
lons_lats_model = np.column_stack((lon.flatten(),lat.flatten()))

maskSP = crmask(sparea,lons_lats_model,lon)
maskST = crmask(starea,lons_lats_model,lon)
areaSP = gr.areacello.where(maskSP==1)
areaST = gr.areacello.where(maskST==1)
areaSP = areaSP.where(landmask==0)
areaST = areaST.where(landmask==0)

sstSPmean = areawgtmean(dfs,gr,maskSP)
sstSTmean = areawgtmean(dfs,gr,maskST)

tempSPmean = areamean(sstSPmean,areaSP)
tempSTmean = areamean(sstSTmean,areaST)

ds = tempSPmean.to_dataset(name='sstSP')
ds.to_netcdf('CESM_tempSP_ssp585.nc')
ds = tempSTmean.to_dataset(name='sstST')
ds.to_netcdf('CESM_tempST_ssp585.nc')

# Salinity
dfs = read_data('gs://cmip6/ScenarioMIP/NCAR/CESM2/ssp585/r1i1p1f1/Omon/so/gn/')

sssSPmean = areawgtmeansalt(dfs,gr,maskSP)
sssSTmean = areawgtmeansalt(dfs,gr,maskST)

saltSPmean = areamean(sssSPmean,areaSP)
saltSTmean = areamean(sssSTmean,areaST)

ds = saltSPmean.to_dataset(name='sssSP')
ds.to_netcdf('CESM_saltSP_ssp585.nc')
ds = saltSTmean.to_dataset(name='sssST')
ds.to_netcdf('CESM_saltST_ssp585.nc')

# MOC
moc = read_data('gs://cmip6/ScenarioMIP/NCAR/CESM2/ssp585/r1i1p1f1/Omon/msftmz/gn/')
amoc = moc.msftmz[:,0,:,:]
amoc26 = amoc.sel(lat=26, method='nearest')
amoc26 = amoc26.groupby('time.year').mean('time')
amoc26time = amoc26.max('lev')*1e-9
amoc26time.to_netcdf('CESM_amoc_ssp585.nc')





# thetao = read_data('gs://cmip6/CMIP/IPSL/IPSL-CM6A-LR/historical/r10i1p1f1/Omon/thetao/gn/')


# df = pd.read_csv('https://storage.googleapis.com/pangeo-cmip6/pangeo-cmip6-zarr-consolidated-stores.csv')
# df.head()

                   #    'activity_id'+'institution_id'+'source_id'+experiment_id'+'memeber_id'+'table_id'+'vari
# able_id'+'grid_label'
# pr = read_data('gs://cmip6/ScenarioMIP/NCAR/CESM2/ssp585/r1i1p1f1/Amon/pr/gn/')

# df_ssp585_pr = df[(df.experiment_id == 'ssp585') & (df.variable_id == 'pr') & (df.institution_id == 'NCAR') &
#  (df.source_id == 'CESM2') & (df.table_id == 'Amon')]


