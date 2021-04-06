import gcsfs

def read_data(uri):
    gcs = gcsfs.GCSFileSystem(token='anon')
    ds = xr.open_zarr(gcs.get_mapper(uri), consolidated=True)
    return ds

##################################################################################

def extract_bs_data(df, year1, year2):
    ds = df.thetao[:,0,130:140,28:43]
    ds = ds.sel(time=slice(year1, year2))
    ds = ds.mean('time')
    ds = ds.to_dataset(name='sst')
    return ds


def extract_bs_data_cnrm(df, year1, year2):
    ds = df.thetao[:,0,212:224,314:330]
    ds = ds.sel(time=slice(year1, year2))
    ds = ds.mean('time')
    ds = ds.to_dataset(name='sst')
    return ds

def extract_bs_data_ec_earth(df, year1, year2):
    ds = df.thetao[:,0,210:222,314:330]
    ds = ds.sel(time=slice(year1, year2))
    ds = ds.mean('time')
    ds = ds.to_dataset(name='sst')
    return ds

###############


# NCAR model
sst_raw_hist = read_data('gs://cmip6/CMIP/NCAR/CESM2-WACCM/historical/r1i1p1f1/Omon/thetao/gr/')
ds = extract_bs_data(sst_raw_hist, '1990', '2010')

sst_raw = read_data('gs://cmip6/ScenarioMIP/NCAR/CESM2-WACCM/ssp585/r1i1p1f1/Omon/thetao/gr/')
ds1 = extract_bs_data(sst_raw, '2040', '2060')
ds2 = extract_bs_data(sst_raw, '2080', '2100')

ds.to_netcdf('CESM2-WACCM_historical_SST_1990_2010.nc')
ds1.to_netcdf('CESM2-WACCM_ssp585_SST_2040_2060.nc')
ds2.to_netcdf('CESM2-WACCM_ssp585_SST_2080_2100.nc')

# GFDL model
sst_raw_hist = read_data('gs://cmip6/CMIP/NOAA-GFDL/GFDL-ESM4/historical/r1i1p1f1/Omon/thetao/gr/')
ds = extract_bs_data(sst_raw_hist, '1990', '2010')

sst_raw = read_data('gs://cmip6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/ssp585/r1i1p1f1/Omon/thetao/gr/')
ds1 = extract_bs_data(sst_raw, '2040', '2060')
ds2 = extract_bs_data(sst_raw, '2080', '2100')
ds.to_netcdf('NOAA-GFDL_historical_SST_1990_2010.nc')
ds1.to_netcdf('NOAA-GFDL_ssp585_SST_2040_2060.nc')
ds2.to_netcdf('NOAA-GFDL_ssp585_SST_2080_2100.nc')

# MRI
sst_raw_hist = read_data('gs://cmip6/CMIP/MRI/MRI-ESM2-0/historical/r1i1p1f1/Omon/thetao/gr/')
ds = extract_bs_data(sst_raw_hist, '1990', '2010')

sst_raw = read_data('gs://cmip6/ScenarioMIP/MRI/MRI-ESM2-0/ssp585/r1i1p1f1/Omon/thetao/gr/')
ds1 = extract_bs_data(sst_raw, '2040', '2060')
ds2 = extract_bs_data(sst_raw, '2080', '2100')
ds.to_netcdf('MRI_historical_SST_1990_2010.nc')
ds1.to_netcdf('MRI_ssp585_SST_2040_2060.nc')
ds2.to_netcdf('MRI_ssp585_SST_2080_2100.nc')

# CNRM
sst_raw_hist = read_data('gs://cmip6/CMIP/CNRM-CERFACS/CNRM-ESM2-1/historical/r1i1p1f2/Omon/thetao/gn/')
ds = extract_bs_data_cnrm(sst_raw_hist, '1990', '2010')

sst_raw = read_data('gs://cmip6/ScenarioMIP/CNRM-CERFACS/CNRM-ESM2-1/ssp585/r1i1p1f2/Omon/thetao/gn/')
ds1 = extract_bs_data_cnrm(sst_raw, '2040', '2060')
ds2 = extract_bs_data_cnrm(sst_raw, '2080', '2100')
ds.to_netcdf('CNRM_historical_SST_1990_2010.nc')
ds1.to_netcdf('CNRM_ssp585_SST_2040_2060.nc')
ds2.to_netcdf('CNRM_ssp585_SST_2080_2100.nc')

# dd=read_data('gs://cmip6/CMIP/NASA-GISS/GISS-E2-1-G-CC/historical/r1i1p1f1/Omon/thetao/gn/')
