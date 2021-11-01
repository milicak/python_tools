import numpy as np
from scipy import interpolate

dfs = xr.open_dataset('/data/opa/bs-mod/upstream_bs-nrt/inicon/S_utss_fullBox-sdn2015m01.nc')
dft = xr.open_dataset('/data/opa/bs-mod/upstream_bs-nrt/inicon/T_utss_fullBox-sdn2015m01.nc')


ds1 = dfs.vosaline[:,:,100,390]
dt1 = dft.votemper[:,:,100,390]

longitude = 37.0*np.ones((3,)) 
latitude = 43.0*np.ones((3,)) 

# cp init_PAPASTATION_m06d15.nc init_BS_NRT.nc
df = xr.open_dataset('init_BS_NRT.nc')

zr_BS = np.copy(dft.deptht)
zr_1D = np.copy(df.deptht)


f = interpolate.interp1d(zr_BS, np.copy(ds1), "nearest", fill_value="extrapolate")  
Snew = f(zr_1D)  
f = interpolate.interp1d(zr_BS, np.copy(dt1), "nearest", fill_value="extrapolate")  
Tnew = f(zr_1D)  


df.votemper[:,:,0,0] = Tnew 
df.votemper[:,:,0,1] = Tnew 
df.votemper[:,:,0,2] = Tnew 
df.votemper[:,:,1,0] = Tnew 
df.votemper[:,:,1,1] = Tnew 
df.votemper[:,:,1,2] = Tnew 
df.votemper[:,:,1,0] = Tnew 
df.votemper[:,:,2,1] = Tnew 
df.votemper[:,:,2,2] = Tnew 
df.votemper[:,:,2,0] = Tnew 


df.vosaline[:,:,0,0] = Snew 
df.vosaline[:,:,0,1] = Snew 
df.vosaline[:,:,0,2] = Snew 
df.vosaline[:,:,1,0] = Snew 
df.vosaline[:,:,1,1] = Snew 
df.vosaline[:,:,1,2] = Snew 
df.vosaline[:,:,1,0] = Snew 
df.vosaline[:,:,2,1] = Snew 
df.vosaline[:,:,2,2] = Snew 
df.vosaline[:,:,2,0] = Snew 


df.to_netcdf('init_BS_NRT_v2_m06d15.nc')

# ds['time_counter']=ds.time_counter+pd.Timedelta("1631 day")
