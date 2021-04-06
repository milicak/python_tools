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


# NCAR model
rsds = read_data('gs://cmip6/CMIP/NCAR/CESM2-WACCM/historical/r1i1p1f1/Amon/rsds/gn/')
rsus = read_data('gs://cmip6/CMIP/NCAR/CESM2-WACCM/historical/r1i1p1f1/Amon/rsus/gn/')
rlds = read_data('gs://cmip6/CMIP/NCAR/CESM2-WACCM/historical/r1i1p1f1/Amon/rlds/gn/')
rlus = read_data('gs://cmip6/CMIP/NCAR/CESM2-WACCM/historical/r1i1p1f1/Amon/rlus/gn/')
hfss = read_data('gs://cmip6/CMIP/NCAR/CESM2-WACCM/historical/r1i1p1f1/Amon/hfss/gn/')
hfls = read_data('gs://cmip6/CMIP/NCAR/CESM2-WACCM/historical/r1i1p1f1/Amon/hfls/gn/')

# qnet=-hfls-hfss+rlds+rsds-rlus-rsus
# https://folk.uib.no/ngfhd/EarthClim/Publications/Papers/Guilyardi_etal_2012.pdf
# A first look at ENSO in CMIP5 Eric Guilyardi
Qnet = rsds.rsds-rsus.rsus+rlds.rlds-rlus.rlus-hfss.hfss-hfls.hfls

Qnet = Qnet.sel(time=slice('1990','2010'))
Qnet = Qnet.mean('time')
Qnet = Qnet.to_dataset(name='Qnet')
Qnet.to_netcdf('CESM2-WACCM_Qnet_1990_2010.nc')

rsds = read_data('gs://cmip6/ScenarioMIP/NCAR/CESM2-WACCM/ssp585/r1i1p1f1/Amon/rsds/gn/')
rsus = read_data('gs://cmip6/ScenarioMIP/NCAR/CESM2-WACCM/ssp585/r1i1p1f1/Amon/rsus/gn/')
rlds = read_data('gs://cmip6/ScenarioMIP/NCAR/CESM2-WACCM/ssp585/r1i1p1f1/Amon/rlds/gn/')
rlus = read_data('gs://cmip6/ScenarioMIP/NCAR/CESM2-WACCM/ssp585/r1i1p1f1/Amon/rlus/gn/')
hfss = read_data('gs://cmip6/ScenarioMIP/NCAR/CESM2-WACCM/ssp585/r1i1p1f1/Amon/hfss/gn/')
hfls = read_data('gs://cmip6/ScenarioMIP/NCAR/CESM2-WACCM/ssp585/r1i1p1f1/Amon/hfls/gn/')
Qnet = rsds.rsds-rsus.rsus+rlds.rlds-rlus.rlus-hfss.hfss-hfls.hfls

Qnet1 = Qnet.sel(time=slice('2040','2060'))
Qnet1 = Qnet1.mean('time')
Qnet1 = Qnet1.to_dataset(name='Qnet')
Qnet1.to_netcdf('CESM2-WACCM_Qnet_2040_2060.nc')

Qnet2 = Qnet.sel(time=slice('2080','2100'))
Qnet2 = Qnet2.mean('time')
Qnet2 = Qnet2.to_dataset(name='Qnet')
Qnet2.to_netcdf('CESM2-WACCM_Qnet_2080_2100.nc')














rsds = read_data('gs://cmip6/ScenarioMIP/MRI/MRI-ESM2-0/ssp585/r1i1p1f1/Amon/rsds/gn/')
rsus = read_data('gs://cmip6/ScenarioMIP/MRI/MRI-ESM2-0/ssp585/r1i1p1f1/Amon/rsus/gn/')
rlds = read_data('gs://cmip6/ScenarioMIP/MRI/MRI-ESM2-0/ssp585/r1i1p1f1/Amon/rlds/gn/')
rlus = read_data('gs://cmip6/ScenarioMIP/MRI/MRI-ESM2-0/ssp585/r1i1p1f1/Amon/rlus/gn/')
hfss = read_data('gs://cmip6/ScenarioMIP/MRI/MRI-ESM2-0/ssp585/r1i1p1f1/Amon/hfss/gn/')
hfls = read_data('gs://cmip6/ScenarioMIP/MRI/MRI-ESM2-0/ssp585/r1i1p1f1/Amon/hfls/gn/')
# Net down surface flux
