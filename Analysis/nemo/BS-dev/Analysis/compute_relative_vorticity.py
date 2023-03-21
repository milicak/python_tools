import numpy as np

gr = xr.open_dataset('/work/opa/da01720/Experiments/CMEMS2_n4.2/bs-test-int-bdy7/model/domain_cfg.nc')

ls1=sorted(glob.glob('/work/opa/da01720/Experiments/CMEMS2_n4.2/bs-test-int-bdy7/rebuilt/*grid_U*'))
ls1=sorted(glob.glob('/work/opa/mi19918/Projects/nemo/BS/BS-NRT_MI_nemo4_2_06/rebuilt/*grid_U*'))
dfu = xr.open_dataset(ls1[-1]) 
ls1=sorted(glob.glob('/work/opa/da01720/Experiments/CMEMS2_n4.2/bs-test-int-bdy7/rebuilt/*grid_V*'))
ls1=sorted(glob.glob('/work/opa/mi19918/Projects/nemo/BS/BS-NRT_MI_nemo4_2_06/rebuilt/*grid_V*'))
dfv = xr.open_dataset(ls1[-1]) 


gr = gr.isel(t=0)


uvel = dfu.uo[0,0,:,:]*gr.e1u
vvel = dfv.vo[0,0,:,:]*gr.e2v

du = uvel[1:,:]-uvel[:-1,:]
du = du.to_dataset(name='du')
tmp = du.du[-1,:]
tmp = tmp.to_dataset(name='du')
du = xr.concat([du,tmp],dim='y')

dv = vvel[:,1:]-vvel[:,:-1]
dv = dv.to_dataset(name='dv')
tmp = dv.dv[:,-1]
tmp = tmp.to_dataset(name='dv')
dv = xr.concat([dv,tmp],dim='x')

zeta = 1/(gr.e1f*gr.e2f) * (dv.dv - du.du)

