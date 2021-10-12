import numpy as np
import xarray as xr


gr = xr.open_dataset('/data/opa/bs-mod/upstream_bs-nrt/geodta/mesh_mask_bs-nrt_s5_smth_changeBosp_Marm_sill_nemo.nc')
df = xr.open_dataset('/work/opa/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/monthly_clim/uTSS_OBC_monthly_v_velocity_nemo_levels_10.nc')
zw = np.copy(gr.gdepw_1d[0,:]) 
dz = np.copy(gr.e3t_0[0,:,0,0])
mask = xr.where(df.v_velocity>-10,1,0)
dnm = mask*dz    
bathy = dnm.sum('z')  
bathy = bathy.to_dataset(name='bathy')    
bathy.to_netcdf('uTSS_interp_bathy.nc') 
