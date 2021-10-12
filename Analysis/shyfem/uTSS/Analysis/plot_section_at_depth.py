import numpy as np
import xarray as xr

def get_data(df, ind_depth):
    elem = np.array(np.copy(df.element_index[0,:,:]-1),dtype=np.int64) 
    ds = df.u_velocity[0,:,ind_depth]  
    level_data = np.copy(ds)   
    level_data[level_data==0]=np.nan 
    d = level_data[elem].mean(axis=1)   
    no_nan_triangles = np.invert(np.isnan(d))   
    elem_no_nan = elem[no_nan_triangles,:] 

    return level_data, elem_no_nan

df = xr.open_dataset('/work/opa/mi19918/Projects/uTSS/Exp_2016_analysis_newTSIC/monthly_clim/uTSS_UV_monthly_clim.nc') 
# Level 66 is depth = 300 meter
level_data, elem_no_nan = get_data(df,66)

plt.tripcolor(df.longitude[0,:],df.latitude[0,:],elem_no_nan[::], level_data,cmap='needJet2');plt.colorbar();
plt.xlim(26,30);plt.ylim(40.5,41.5); 


plt.tripcolor(df.longitude[0,:],df.latitude[0,:],elem_no_nan3[::], level_data3,cmap='needJet2',edgecolors='k');plt.color
          ...: bar();plt.xlim(27,29.5);plt.ylim(40.5,41.); plt.clim(-1e-2,1e-2)
