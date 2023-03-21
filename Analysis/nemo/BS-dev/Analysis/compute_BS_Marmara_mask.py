import numpy as np
from matplotlib.path import Path

gr = xr.open_dataset('/work/opa/da01720/Experiments/CMEMS2_n4.2/bs-test-int-bdy7/model/domain_cfg.nc')

lonrec = np.array([28.25,29.3,29.3,28.25,28.25])
latrec = np.array([40.68,40.68,41.1,41.1,40.68])
vertices = np.vstack((lonrec,latrec))
mpath = Path( np.transpose(vertices) )

lonlat = np.dstack((np.copy(gr.nav_lon), np.copy(gr.nav_lat))) 
lonlat_flat = lonlat.reshape((-1, 2))   

mask_flat = mpath.contains_points(lonlat_flat)     
mask = mask_flat.reshape(gr.nav_lon.shape)       

ds = xr.Dataset(
    data_vars=dict(
        mask=(["y", "x"], mask),
    ),
)

ds['mask'] = np.abs(ds.mask-1)
ds.to_netcdf('BlackSea_mask.nc')
