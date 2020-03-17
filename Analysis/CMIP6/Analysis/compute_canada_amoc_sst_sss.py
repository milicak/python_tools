import gcsfs
from matplotlib import path

def read_data(uri):
    gcs = gcsfs.GCSFileSystem(token='anon')
    ds = xr.open_zarr(gcs.get_mapper(uri), consolidated=True)
    return ds


def areamean(var,area, lonvar='nlat', latvar='nlon'):
    ds = var.sum((latvar,lonvar),skipna=True)/area.sum((latvar,lonvar),skipna=True)
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
    sssmean = sssmean*gr.areacello
    return sssmean



# Subpolar gyre area
sparea = np.genfromtxt('subpolar_gyre_lonlat.txt',usecols=(0,1))
# Subtropical gyre area
starea = np.genfromtxt('subtropical_gyre_lonlat.txt',usecols=(0,1))

dfs = read_data('gs://cmip6/ScenarioMIP/CCCma/CanESM5/ssp585/r1i1p1f1/Omon/tos/gn/')
gr = xr.open_dataset('areacello_Ofx_CanESM5_ssp585_r1i1p1f1_gn.nc')

lon = np.copy(dfs.longitude)
lat = np.copy(dfs.latitude)
lons_lats_model = np.column_stack((lon.flatten(),lat.flatten()))

maskSP = crmask(sparea,lons_lats_model,lon)
maskST = crmask(starea,lons_lats_model,lon)

areaSP = gr.areacello.where(maskSP==1)
areaST = gr.areacello.where(maskST==1)

sstSPmean = areawgtmean(dfs,gr,maskSP)
sstSTmean = areawgtmean(dfs,gr,maskST)

tempSPmean = areamean(sstSPmean,areaSP,lonvar='i',latvar='j')
tempSTmean = areamean(sstSTmean,areaST,lonvar='i',latvar='j')

tempSPmean.to_netcdf('CanESM_tempSP_ssp585.nc')
tempSTmean.to_netcdf('CanESM_tempST_ssp585.nc')

# Salinity
dfs = read_data('gs://cmip6/ScenarioMIP/CCCma/CanESM5/ssp585/r1i1p1f1/Omon/so/gn/')

sssSPmean = areawgtmeansalt(dfs,gr,maskSP)
sssSTmean = areawgtmeansalt(dfs,gr,maskST)
saltSPmean = areamean(sssSPmean,areaSP,lonvar='i',latvar='j')
saltSTmean = areamean(sssSTmean,areaST,lonvar='i',latvar='j')

saltSPmean.to_netcdf('CanESM_saltSP_ssp585.nc')
saltSTmean.to_netcdf('CanESM_saltST_ssp585.nc')

# MOC
moc = read_data('gs://cmip6/ScenarioMIP/CCCma/CanESM5/ssp585/r1i1p1f1/Omon/msftmz/gn/')
amoc = moc.msftmz[:,0,:,:]
amoc26 = amoc.sel(lat=26, method='nearest')
amoc26 = amoc26.groupby('time.year').mean('time')
amoc26time = amoc26.max('lev')*1e-9
amoc26time.to_netcdf('CanESM_amoc_ssp585.nc')





# thetao = read_data('gs://cmip6/CMIP/IPSL/IPSL-CM6A-LR/historical/r10i1p1f1/Omon/thetao/gn/')


# df = pd.read_csv('https://storage.googleapis.com/pangeo-cmip6/pangeo-cmip6-zarr-consolidated-stores.csv')
# df.head()

                   #    'activity_id'+'institution_id'+'source_id'+experiment_id'+'memeber_id'+'table_id'+'vari
# able_id'+'grid_label'
# pr = read_data('gs://cmip6/ScenarioMIP/NCAR/CESM2/ssp585/r1i1p1f1/Amon/pr/gn/')

# df_ssp585_pr = df[(df.experiment_id == 'ssp585') & (df.variable_id == 'pr') & (df.institution_id == 'NCAR') &
#  (df.source_id == 'CESM2') & (df.table_id == 'Amon')]


