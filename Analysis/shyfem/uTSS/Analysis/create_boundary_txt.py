import numpy as np
import xarray as xr


gr = xr.open_dataset('/work/mi19918/Projects/uTSS/preproc/data/grid/Marmara_basbathy_ser_grd.nc')

boundary_nodes_path = 'bound3.txt'
boundary_nodes_outpath = 'uTSS_boundary_list_lonlat3.txt'

data = np.genfromtxt(boundary_nodes_path,usecols=(0))
aa = np.zeros(data.shape)

for ii, ind in enumerate(data):
    aa[ii] = np.copy(np.where(gr.pointIndexes==ind))


aa = np.int32(aa)
lonlat = np.copy(gr.pointsnew[aa,:])
np.savetxt(boundary_nodes_outpath, lonlat, fmt='%4.8f')
