import xarray as xr
import numpy as np
import xgcm


gr = xr.open_dataset('/work/opa/sc02915/data/bs-simu_6.6_conf/mesh_mask_Aug19_31-5-25.nc')
gr2 = xr.open_dataset('/work/opa/sc02915/data/bs-simu_6.6_conf/bathy_meter_Aug19.nc')   
gr = gr.squeeze('t')

# gr = gr.rename({'t': 'time', 'z': 'depth', 'x': 'lon', 'y': 'lat'})
# glamt      (y_c, x_c) float64 ...
# glamf      (y_f, x_f) float64 ...
# gphit      (y_c, x_c) float64 ...
# gphif      (y_f, x_f) float64 ...
# e1t        (y_c, x_c) float64 ...
# e1u        (y_c, x_f) float64 ...
# e1v        (y_f, x_c) float64 ...
# e1f        (y_f, x_f) float64 ...
# e2t        (y_c, x_c) float64 ...
# e2u        (y_c, x_f) float64 ...
# e2v        (y_f, x_c) float64 ...
# e2f        (y_f, x_f) float64 ...
# e3t_0      (z_c, y_c, x_c) float64 ...
# e3u_0      (z_c, y_c, x_f) float64 ...
# e3v_0      (z_c, y_f, x_c) float64 ...
# e3f_0      (z_c, y_f, x_f) float64 ...
# e3w_0      (z_f, y_c, x_c) float64 ...
# ht_0       (y_c, x_c) float64 ...
# fmaskutil  (y_f, x_f) float64 ...
domcfg = xr.Dataset({'glamt': (('y_c','x_c'),gr.glamt),
                'glamf': (('y_f','x_f'),gr.glamf),
                'gphit': (('y_c','x_c'),gr.gphit),
                'gphif': (('y_f','x_f'),gr.gphif),
                'e1t'  : (('y_c','x_c'),gr.e1t),
                'e1u'  : (('y_c','x_f'),gr.e1u),
                'e1v'  : (('y_f','x_c'),gr.e1v),
                'e1f'  : (('y_f','x_f'),gr.e1f),
                'e2t'  : (('y_c','x_c'),gr.e2t),
                'e2u'  : (('y_c','x_f'),gr.e2u),
                'e2v'  : (('y_f','x_c'),gr.e2v),
                'e2f'  : (('y_f','x_f'),gr.e2f),
                'e3t_0': (('z_c','y_c','x_c'),gr.e3t_0),
                'e3u_0': (('z_c','y_c','x_f'),gr.e3u_0),
                'e3v_0': (('z_c','y_f','x_c'),gr.e3v_0),
                'e3f_0': (('z_c','y_f','x_f'),gr.e3t_0),
                'e3w_0': (('z_f','y_c','x_c'),gr.e3w_0),
                 'ht_0' : (('y_c','x_c'),gr2.Bathymetry[0,:,:]),
                'fmaskutil': (('y_f','x_f'),gr.fmaskutil)},
                    {'x_c': np.copy(gr.glamt[0,:]), 'y_c':
                     np.copy(gr.gphit[:,0]), 
                 'x_f': np.copy(gr.glamf[0,:]), 'y_f': np.copy(gr.gphif[:,0]),
                 'z_c': np.copy(gr.e3t_1d), 'z_f': np.copy(gr.e3t_1d)})

xc=np.arange(0,395,1)  
yc=np.arange(0,215,1)  
xf=xc+0.5
yf=yc+0.5 
domcfg = xr.Dataset({'glamt': (('y_c','x_c'),gr.glamt),
                'glamf': (('y_f','x_f'),gr.glamf),
                'gphit': (('y_c','x_c'),gr.gphit),
                'gphif': (('y_f','x_f'),gr.gphif),
                'e1t'  : (('y_c','x_c'),gr.e1t),
                'e1u'  : (('y_c','x_f'),gr.e1u),
                'e1v'  : (('y_f','x_c'),gr.e1v),
                'e1f'  : (('y_f','x_f'),gr.e1f),
                'e2t'  : (('y_c','x_c'),gr.e2t),
                'e2u'  : (('y_c','x_f'),gr.e2u),
                'e2v'  : (('y_f','x_c'),gr.e2v),
                'e2f'  : (('y_f','x_f'),gr.e2f),
                'e3t_0': (('z_c','y_c','x_c'),gr.e3t_0),
                'e3u_0': (('z_c','y_c','x_f'),gr.e3u_0),
                'e3v_0': (('z_c','y_f','x_c'),gr.e3v_0),
                'e3f_0': (('z_c','y_f','x_f'),gr.e3t_0),
                'e3w_0': (('z_f','y_c','x_c'),gr.e3w_0),
                 'ht_0' : (('y_c','x_c'),gr2.Bathymetry[0,:,:]),
                'fmaskutil': (('y_f','x_f'),np.float32(gr.fmaskutil))},
                {'x_c': np.copy(xc), 'y_c': np.copy(yc), 
                 'x_f': np.copy(xf), 'y_f': np.copy(yf),
                 'z_c': np.copy(gr.e3t_1d), 'z_f': np.copy(gr.e3t_1d)-0.5})

domcfg.x_c.attrs = {'axis': 'X'}
domcfg.y_c.attrs = {'axis': 'Y'}
domcfg.z_c.attrs = {'axis': 'Z'}

# domcfg2.attrs['DOMAIN_number_total']=np.int64(1)
# domcfg2.attrs['DOMAIN_number']=np.int64(1)
# domcfg2.attrs['DOMAIN_dimensions_ids']=np.array([1,2],dtype=np.int32)
# domcfg2.attrs['DOMAIN_size_global']=np.array([20,40],dtype=np.int32)
# domcfg2.attrs['DOMAIN_size_local']=np.array([20,40],dtype=np.int32)
# domcfg2.attrs['DOMAIN_position_first']=np.array([1,1],dtype=np.int32)
# domcfg2.attrs['DOMAIN_position_last']=np.array([20,40],dtype=np.int32)
# domcfg2.attrs['DOMAIN_halo_size_start']=np.array([0,0],dtype=np.int32)
# domcfg2.attrs['DOMAIN_halo_size_end']=np.array([0,0],dtype=np.int32)
# domcfg2.attrs['DOMAIN_type']='BOX'
# domcfg2.attrs['nn_cfg']=np.int32(2)
# domcfg2.attrs['cn_cfg']='BASIN'
#

df=xr.open_dataset('/data/inputs/metocean/historical/model/ocean/CMCC/CMEMS/reanalysis/day/2019/12/20191231_d-CMCC--TEMP-BSe3r1-BS-b20200901_re-sv09.00.nc')


# xgcm properties
bd={'boundary':'fill', 'fill_value':0}

metrics = {
    ('X',): ['e1t', 'e1u', 'e1v', 'e1f'], # X distances
    ('Y',): ['e2t', 'e2u', 'e2v', 'e2f'], # Y distances
    ('Z',): ['e3t_0', 'e3u_0', 'e3v_0', 'e3f_0', 'e3w_0'], # Z distances
    #('X', 'Y'): [] # Areas TODO
}
grid = xgcm.Grid(domcfg, metrics=metrics, periodic=False)
print(grid)




# create dataset
ds = xr.Dataset({'vol_sigma_tr': (('time','lat','sigma2_bin'),vol_sigma_tr)},
                {'e3t': (('z_c', 'y_c', 'x_c'), gr.e3t_0[0,:,:,:])}
                {'time': time, 'lat': df.lat, 'sigma2_bin': df.sigma2_bin})

ds = xr.Dataset({'vol_sigma_tr': (('time','lat','sigma2_bin'),vol_sigma_tr)},
                {'time': time, 'lat': df.lat, 'sigma2_bin': df.sigma2_bin})
