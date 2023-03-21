import xarray as xr

df = xr.open_dataset('/work/opa/mi19918/Projects/uTSS_SHYFEM/work/tmp3/REG_1448_1449_0_salinity_0.0005.nc')
df=xr.open_dataset('/work/opa/mi19918/Projects/uTSS_SHYFEM/work/tmp3/REG_1558_1559_0_salinity_0.0005.nc')
df=xr.open_dataset('/work/opa/mi19918/Projects/uTSS_SHYFEM/work/tmp3/REG_1669_1670_0_salinity_0.0005.nc')
df.salinity[0,:,180:290,60:170] = np.nan
# get the surface value
d1 = df.salinity[0,0,:,:]

# get the lonmean values
lonmean = np.zeros(df.lat.shape[0])
saltmean = np.zeros((df.level.shape[0],df.lat.shape[0]))
for ind in range(0,lonmean.shape[0]):
    # logic
    dnm = d1[ind,:].notnull()
    lonmean[ind] = d1[ind,dnm].lon.mean()
    saltmean[:,ind] = df.salinity[0,:,ind,:].interp(lon=lonmean[ind], method="nearest")


# change southern part
lonmean[0:160] = 29.0
