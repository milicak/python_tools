import numpy as np
import xesmf as xe

# Dimensions:       (time_counter: 12, x: 182, y: 149, z: 31)
# Coordinates:
#   * time_counter  (time_counter) float32 1.296e+06 3.888e+06 ... 2.981e+07
# Dimensions without coordinates: x, y, z
# Data variables:
#     nav_lon       (y, x) float32 ...
#     nav_lat       (y, x) float32 ...
#     depth         (z) float32 ...
#     vosaline      (time_counter, z, y, x) float32 ...


# lets use first initial conditions 
times = pd.date_range("2000-01-01", "2000-12-31", name="time_counter",freq='1M')

df1 = xr.open_dataset('/data/opa/bs-mod/upstream_bs-nrt/baseline/inicon/S_utss_fullBox-sdn2015m01.nc')
df2 = xr.open_dataset('/data/opa/bs-mod/upstream_bs-nrt/baseline/inicon/T_utss_fullBox-sdn2015m01.nc')
df1['time_counter'] = np.array([times[0]])
df2['time_counter'] = np.array([times[0]])

ds1 = df1.reindex(time_counter=times, method='ffill')
ds2 = df2.reindex(time_counter=times, method='ffill')

ls1 = sorted(glob.glob('/work/opa/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/monthly_clim/*salinity*nemo_level*'))
dd1 = xr.open_mfdataset(ls1,concat_dim='time', combine='nested')
ls1 = sorted(glob.glob('/work/opa/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/monthly_clim/*temperature*nemo_level*'))
dd2 = xr.open_mfdataset(ls1,concat_dim='time', combine='nested')

xind1 = 40
xind2 = 85
yind1 = 5
yind2 = 22

# nemo grid
gr1 = df1['lon'][0,:]          
gr1 = gr1.to_dataset(name='lon')   
gr1['lat']=df1['lat'][:,0]
# shyfem grid
gr2 = dd['lon'][0,0,:].load()          
gr2 = gr2.to_dataset(name='lon')   
gr2['lat']=dd['lat'][0,:,0].load()

regridder = xe.Regridder(gr2, gr1, 'bilinear') 
#apply regridder                  
dr_out1 = regridder(dd1['salinity'])  
dr_out2 = regridder(dd2['temperature'])  
ds1.vosaline[:,:,yind1:yind2,xind1:xind2] = dr_out1[:,:,yind1:yind2,xind1:xind2]
ds2.votemper[:,:,yind1:yind2,xind1:xind2] = dr_out2[:,:,yind1:yind2,xind1:xind2]

ds1.to_netcdf('/work/opa/mi19918/Projects/shared/marmara_salt_nudging.nc',unlimited_dims='time_counter')
ds2.to_netcdf('/work/opa/mi19918/Projects/shared/marmara_temp_nudging.nc',unlimited_dims='time_counter')

