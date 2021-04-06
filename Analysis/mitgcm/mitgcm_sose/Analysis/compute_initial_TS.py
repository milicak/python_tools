import numpy as np
from scipy.io import loadmat

dfsose = loadmat('/okyanus/users/milicak/dataset/SOSE/sosemonthly_uvts_nan.mat')
saltc = dfsose['saltc']
tempc = dfsose['tempc']

grsose = loadmat('/okyanus/users/milicak/Analysis/mitgcm/mitgcm_sose/Analysis/sosegrid_mitgcm.mat')
lonc = grsose['XC']
latc = grsose['YC']
depth = grsose['RC']


ds = xr.DataArray(saltc,dims=('lon', 'lat', 'depth'),
                  coords={'lon': lonc[:,0], 'lat': latc[0,:], 'depth':
                          depth[:,0]})
ds_salt = ds.to_dataset(name='saltc')
ds = xr.DataArray(tempc,dims=('lon', 'lat', 'depth'),
                  coords={'lon': lonc[:,0], 'lat': latc[0,:], 'depth':
                          depth[:,0]})
ds_temp = ds.to_dataset(name='tempc')

# x = xr.DataArray(lonc,dims={'lon', 'lat'})
# y = xr.DataArray(latc,dims={'lon', 'lat'})
# ds_temp2 = ds_temp.interp(lon=x[:,0],lat=y[0,:])
# ds_salt2 = ds_salt.interp(lon=x[:,0],lat=y[0,:])

x = xr.open_dataset('newgrid_XC_info.nc')
y = xr.open_dataset('newgrid_YC_info.nc')
ds_temp2 = ds_temp.interp(lon=x.XC,lat=y.YC)
ds_salt2 = ds_salt.interp(lon=x.XC,lat=y.YC)




