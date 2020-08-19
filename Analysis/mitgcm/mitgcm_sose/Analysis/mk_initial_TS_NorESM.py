import numpy as np
import xesmf as xe
import xmitgcm


df = xmitgcm.open_mdsdataset('/archive/milicak/MITgcm_c65/Projects/mitgcm_sose/ctrl/',
                                prefix='THETA'
                                ,geometry='sphericalpolar',
                                ref_date='2003-12-31 12:0:0',delta_t=80)


sst_model = df.THETA.isel(Z=0)
sst_model = sst_model.rename({'XC':'lon', 'YC':'lat'})
sst_model = sst_model.mean('time').load()

# lon_mitgcm = np.copy(sst_model['lon'])
# lat_mitgcm = np.copy(sst_model['lat'])
# lon_mitgcm,lat_mitgcm=np.meshgrid(lon_mitgcm,lat_mitgcm)
# lon_mitgcm=np.transpose(lon_mitgcm)
# lat_mitgcm=np.transpose(lat_mitgcm)
# mitgcm_grid = xr.merge([lon_mitgcm,lat_mitgcm])
#
# NorESM grid
grid = xr.open_dataset('/archive/milicak/dataset/NorESM/grid.nc')
gr = grid.rename({'x':'lon', 'y': 'lat'})
#
lat_noresm = gr.plat
lon_noresm = gr.plon

NorESM_grid = xr.merge([lon_noresm,lat_noresm])


NorESM_grid=NorESM_grid.rename_vars({'plat':'lat','plon':'lon'})

regridder = xe.Regridder(NorESM_grid, sst_model, 'bilinear', periodic=True)


data = xr.open_dataset('/archive/milicak/MITgcm_c65/Projects/mitgcm_sose/ctrl/NOIIAJRAOC20TR_TL319_tn14_20190709.micom.hm.1663-01.nc')

data_out = regridder2(data['sst'])

