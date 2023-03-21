import numpy as np                                                              
import xarray as xr                                                             
                                                                                
df = xr.open_dataset('/work/opa/mi19918/Projects/nemo/BS/BS-NRT_MI_bdy7_geom3_rlxbiharm02//model/eddy_viscosity_2D.nc')
gr = xr.open_dataset('/work/opa/mi19918/Projects/nemo/BS/BS-NRT_MI_bdy7_geom3/model/domain_cfg.nc')

# 1/12 Uv*Lv^3 
rn_Uv = 0.05
rn_Lv = 10e2
cff1 = 0.07/rn_Uv

# biharm = (rn_Uv*rn_Lv**3)/12
biharm = (cff1*rn_Uv*gr.e2t**3)/12

# increase tip of danube!
biharm[0,180:190,95:104] = 10*biharm[0,180:190,95:104]
# increase bosphorus strait
# biharm[0,19:30,66:78] = 2*biharm[0,19:30,66:78]


biharm = biharm.rename({'t': 'TIME_COUNTER', 'y': 'Y', 'x': 'X'})
biharm = biharm.where(df.ahmt_2d!=0)
biharm = biharm.fillna(0)
# make it 3D
tmp = biharm.expand_dims(dim={"z": 121},axis=1)
tmp = tmp.rename({'TIME_COUNTER': 't', 'Y': 'y', 'X': 'x'})
df['ahmt_2d'] = biharm
df['ahmf_2d'] = biharm
df.to_netcdf('eddy_viscosity_2D.nc')                                              
df.close()
               
# 3D version
variables = ['e3t_0','e3f_0']
ds = gr[variables]
ds = ds.rename({'e3t_0':'ahmt_3d','e3f_0':'ahmf_3d'})
ds['ahmt_3d'] = tmp
ds['ahmf_3d'] = tmp
