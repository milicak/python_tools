import gcsfs

def read_data(uri):
    gcs = gcsfs.GCSFileSystem(token='anon')
    ds = xr.open_zarr(gcs.get_mapper(uri), consolidated=True)
    return ds

##################################################################################

def extract_bs_data(df, area, year1, year2):
    ds = df.thetao[:,0,130:140,28:43]
    ds = ds.groupby('time.year').mean('time')
    ds = ds.sel(year=slice(year1, year2))
    ds = ds * area[130:140,28:43]
    tmparea = area[130:140,28:43]
    ds = ds.sum(['lon','lat'])/tmparea.sum(['lon','lat'])
    ds = ds.to_dataset(name='sst')
    return ds


def extract_bs_data_cnrm(df, area, year1, year2):
    ds = df.thetao[:,0,212:224,314:330]
    ds = ds.groupby('time.year').mean('time')
    ds = ds.sel(year=slice(year1, year2))
    ds = ds * area[212:224,314:330]
    tmparea = area[212:224,314:330]
    ds = ds.sum(['x','y'])/tmparea.sum(['x','y'])
    ds = ds.to_dataset(name='sst')
    return ds

###############

gra = xr.open_dataset('/okyanus/users/dcetin/BS/data/areacello_Ofx_GFDL-CM4_piControl_r1i1p1f1_gr.nc')
area = gra.areacello
gr = xr.open_dataset('areacello_Ofx_CNRM-CM6-1_historical_r11i1p1f2_gn.nc')
area2 = gr.areacello

# NCAR model
sst_raw_hist = read_data('gs://cmip6/CMIP/NCAR/CESM2-WACCM/historical/r1i1p1f1/Omon/thetao/gr/')
ds = extract_bs_data(sst_raw_hist, area, '1900', '2014')
ds.to_netcdf('CESM2-WACCM_historical_SST.nc')

sst_raw = read_data('gs://cmip6/ScenarioMIP/NCAR/CESM2-WACCM/ssp585/r1i1p1f1/Omon/thetao/gr/')
ds = extract_bs_data(sst_raw, area, '2015', '2100')
ds.to_netcdf('CESM2-WACCM_ssp585_SST.nc')

# GFDL model
sst_raw_hist = read_data('gs://cmip6/CMIP/NOAA-GFDL/GFDL-ESM4/historical/r1i1p1f1/Omon/thetao/gr/')
ds = extract_bs_data(sst_raw_hist, area, '1900', '2014')
ds.to_netcdf('NOAA-GFDL_historical_SST.nc')

sst_raw = read_data('gs://cmip6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp585/r1i1p1f1/Omon/thetao/gr/')
ds = extract_bs_data(sst_raw, area, '2015', '2100')
ds.to_netcdf('NOAA-GFDL_ssp585_SST.nc')

# MRI
sst_raw_hist = read_data('gs://cmip6/CMIP/MRI/MRI-ESM2-0/historical/r1i1p1f1/Omon/thetao/gr/')
ds = extract_bs_data(sst_raw_hist, area, '1900', '2014')
ds.to_netcdf('MRI_historical_SST.nc')

sst_raw = read_data('gs://cmip6/ScenarioMIP/MRI/MRI-ESM2-0/ssp585/r1i1p1f1/Omon/thetao/gr/')
ds = extract_bs_data(sst_raw, area, '2015', '2100')
ds.to_netcdf('MRI_ssp585_SST.nc')

# CNRM
sst_raw_hist = read_data('gs://cmip6/CMIP/CNRM-CERFACS/CNRM-ESM2-1/historical/r1i1p1f2/Omon/thetao/gn/')
ds = extract_bs_data_cnrm(sst_raw_hist, area2, '1900', '2014')
ds.to_netcdf('CNRM_historical_SST.nc')

sst_raw = read_data('gs://cmip6/ScenarioMIP/CNRM-CERFACS/CNRM-ESM2-1/ssp585/r1i1p1f2/Omon/thetao/gn/')
ds = extract_bs_data_cnrm(sst_raw, area2, '2015', '2100')
ds.to_netcdf('CNRM_ssp585_SST.nc')


