import xclim
import xclim.utils
from xclim import icclim
import xarray as xr
from dask.distributed import Client

client = Client(n_workers=10)
client

root_folder = '/archive/milicak/dataset/CMIP6_nas/Raw_Data/day/'
root_folder2 = '/archive/milicak/dataset/CMIP6/indices/'

expid = 'GFDL-ESM4'

fname1 = root_folder + 'pr/' + expid + '/pr_day_GFDL-ESM4_hist_ssp585_gr1_185001-210012.nc'
pr = xr.open_dataset(fname1)
# pr = xr.open_dataset(fname1, chunks={'time': 10000})
fname2 = root_folder + 'tasmax/' + expid + '/tasmax_day_GFDL-ESM4_hist_ssp585_gr1_185001-210012.nc'
tmax = xr.open_dataset(fname2)
# tmax = xr.open_dataset(fname2, chunks={'time': 10000})
fname3 = root_folder + 'tasmin/' + expid + '/tasmin_day_GFDL-ESM4_hist_ssp585_gr1_185001-210012.nc'
tmin = xr.open_dataset(fname3)
# tmin = xr.open_dataset(fname3, chunks={'time': 10000})
fname4 = root_folder + 'tas/' + expid + '/tas_day_GFDL-ESM4_hist_ssp585_gr1_185001-210012.nc'
tas = xr.open_dataset(fname4)
# tas = xr.open_dataset(fname4, chunks={'time': 10000})


# Compute indices
frost_days = xclim.indices.frost_days(tmin.tasmin, freq="YS")
frost_days = frost_days.to_dataset(name='frost')
outname = root_folder2 + expid + '/frost_days.nc'
frost_days.to_netcdf(outname)

summer_days = xclim.indices.tx_days_above(tmax.tasmax, freq="YS")
summer_days = summer_days.to_dataset(name='summer_days')
outname = root_folder2 + expid + '/summer_days.nc'
summer_days.to_netcdf(outname)


icing_days = xclim.indices.ice_days(tmax.tasmax, freq="YS")
icing_days = icing_days.to_dataset(name='icing_days')
outname = root_folder2 + expid + '/icing_days.nc'
icing_days.to_netcdf(outname)

tropical_days = xclim.indices.tropical_nights(tmin.tasmin, freq="YS")
tropical_days = tropical_days.to_dataset(name='tropical_days')
outname = root_folder2 + expid + '/tropical_days.nc'
tropical_days.to_netcdf(outname)

growing_season_length = xclim.indices.growing_season_length(tas.tas, freq="YS")
growing_season_length = growing_season_length.to_dataset(name='growing_season_length')
outname = root_folder2 + expid + '/growing_season_length.nc'
growing_season_length.to_netcdf(outname)

# Monthly maximum value of daily maximum temperature TXx
txx = xclim.indices.tx_max(tmax.tasmax, freq="YS")
txx = txx.to_dataset(name='txx')
outname = root_folder2 + expid + '/txx.nc'
txx.to_netcdf(outname)

# Monthly maximum value of daily minimum temperature TNx
tnx = xclim.indices.tn_max(tmin.tasmin, freq="YS")
tnx = tnx.to_dataset(name='tnx')
outname = root_folder2 + expid + '/tnx.nc'
tnx.to_netcdf(outname)

# Monthly minimum value of daily maximum temperature TXn
txn = xclim.indices.tx_min(tmax.tasmax, freq="YS")
txn = txn.to_dataset(name='txn')
outname = root_folder2 + expid + '/txn.nc'
txn.to_netcdf(outname)

# Monthly minimum value of daily minimum temperature TNn
tnn = xclim.indices.tn_min(tmin.tasmin, freq="YS")
tnn = tnn.to_dataset(name='tnn')
outname = root_folder2 + expid + '/tnn.nc'
tnn.to_netcdf(outname)

df = tmin.sel(time=slice('1980','2100'))
# t10 = xclim.utils.percentile_doy(tmin.tasmin, per=0.1)
# tn10p = icclim.TN10p(tmin.tasmin, t10, freq = "YS")
t10 = xclim.utils.percentile_doy(df.tasmin, per=0.1)
tn10p = icclim.TN10p(df.tasmin, t10, freq = "YS")
tn10p = tn10p.to_dataset(name='tn10p')
outname = root_folder2 + expid + '/tn10p.nc'
tn10p.to_netcdf(outname)

t90 = xclim.utils.percentile_doy(df.tasmin, per=0.9)
tn90p = icclim.TN90p(df.tasmin, t90, freq = "YS")
tn90p = tn90p.to_dataset(name='tn90p')
outname = root_folder2 + expid + '/tn90p.nc'
tn90p.to_netcdf(outname)

df = tmax.sel(time=slice('1980','2100'))
t10 = xclim.utils.percentile_doy(df.tasmax, per=0.1)
tx10p = icclim.TX10p(df.tasmax, t10, freq = "YS")
tx10p = tx10p.to_dataset(name='tx10p')
outname = root_folder2 + expid + '/tx10p.nc'
tx10p.to_netcdf(outname)

t90 = xclim.utils.percentile_doy(df.tasmax, per=0.9)
tx90p = icclim.TX90p(df.tasmax, t90, freq = "YS")
tx90p = tx90p.to_dataset(name='tx90p')
outname = root_folder2 + expid + '/tx90p.nc'
tx90p.to_netcdf(outname)

wsdi = xclim.indices.warm_spell_duration_index(tmax.tasmax, t90, freq = "YS")
wsdi = wsdi.to_dataset(name='wsdi')
outname = root_folder2 + expid + '/wsdi.nc'
wsdi.to_netcdf(outname)



csdi = xclim.indices.cold_spell_duration_index(tmin.tasmin, freq = "YS")
dtr = xclim.indices.daily_temperature_range(tmax.tasmax, tmin.tasmin, freq = "YS")
rx1day = xclim.indices.max_1day_precipitation_amount(pr.pr, freq = "YS")
window = 5 # 5 day window
rx5day = xclim.max_n_day_precipitation_amount(pr.pr, window=window, freq="YS")
sdii = xclim.indices.daily_pr_intensity(pr.pr, thresh='1 mm/day', freq="YS")

r10mm = xclim.icclim.R10mm(pr.pr, freq="YS")
r20mm = xclim.icclim.R20mm(pr.pr, freq="YS")
cdd = xclim.indices.maximum_consecutive_dry_days(pr.pr, thresh='1 mm/day', freq="YS")
cwd = xclim.indices.maximum_consecutive_wet_days(pr.pr, thresh='1 mm/day', freq="YS")
# r95p =
# r99p =
# prcptot


