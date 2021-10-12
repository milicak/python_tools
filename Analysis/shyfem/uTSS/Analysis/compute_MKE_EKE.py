import numpy as np
import xrscipy.signal as dsp  



ls1 = sorted(glob.glob('/work/opa/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/OUT/*ous*'))
df = xr.open_mfdataset(ls1)
dfu = df.u_velocity[:,:,0]
dfv = df.v_velocity[:,:,0]
ndays = dfu.shape[0]
t = np.linspace(1, ndays, ndays)
dfu['time'] = t
dfv['time'] = t

# low pass filter frequency fr = 1/days i.e. 50 days filter is 1/50
fr = 1.0/50.0
sig = xr.DataArray(np.sin(0.016*t) + np.random.normal(0, 0.1, t.size),coords=[('time', t)], name='signal') 
lowu = dsp.lowpass(dfu.load(), fr, order=8, dim='time')
lowv = dsp.lowpass(dfv.load(), fr, order=8, dim='time')

upr = dfu - lowu
vpr = dfv - lowv

MKE = 0.5*(lowu**2+lowv**2)
EKE = 0.5*(upr**2+vpr**2)

MKE = MKE.to_dataset(name='MKE')
EKE = EKE.to_dataset(name='EKE')

MKE.to_netcdf('/work/opa/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/daily_clim/MKE_daily.nc')
EKE.to_netcdf('/work/opa/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/daily_clim/EKE_daily.nc')
