import numpy as np


root_folder = '/work/opa/mi19918/Projects/nemo/BS/MOC_data/'
filenames = root_folder + 'sigma2_BSe3r1*'

ls1 = sorted(glob.glob(filenames))
sigma2 = np.zeros((31,215,395))
for ind in range(9861,len(ls1)):
    print(ind)
    df = xr.open_dataset(ls1[ind])
    sigma2 = sigma2 + np.copy(df.sigma2)     
    if np.nansum(sigma2)==0:
        print('mehmet',ind)

sigma2 = sigma2/366 #len(ls1)
ds = xr.open_dataset('/work/opa/mi19918/Projects/nemo/BS/MOC_data/sigma2_mean_time.nc')
sigma2 = 0.5*(sigma2+np.copy(ds.sigma2))
dd = xr.Dataset({'sigma2': (('depth','lat','lon'),sigma2)},
                {'depth': ds.depth, 'lat': ds.lat, 'lon': ds.lon})

dd.to_netcdf('/work/opa/mi19918/Projects/nemo/BS/MOC_data/sigma2_mean_time_1993_2020.nc')




ls1 = sorted(glob.glob('/work/opa/mi19918/Projects/nemo/BS/MOC_data/*meridional*sigma2*BSe3r1*')) 
ds = xr.open_dataset('/users_home/opa/mi19918/moc_meridional_sigma2_BSe3r1_mean.nc')
vol_sigma_tr = np.zeros((215,100))
for ind in range(0,len(ls1)):
    print(ind)
    df = xr.open_dataset(ls1[ind])
    vol_sigma_tr = vol_sigma_tr + np.copy(df.vol_sigma_tr)     


vol_sigma_tr = vol_sigma_tr/len(ls1)
dd = xr.Dataset({'vol_sigma_tr': (('lat','sigma2_bin'),-vol_sigma_tr)},
                {'lat': ds.lat, 'sigma2_bin': ds.sigma2_bin})
dd.to_netcdf('/users_home/opa/mi19918/moc_meridional_sigma2_BSe3r1_mean_1993_2020.nc')






# df = xr.open_mfdataset(ls1,concat_dim='time') 
# time = pd.date_range("1993-01-01", freq="D", periods=10227)  
# df['time'] = time
# ds = df.mean('time')
# ds = xr.open_dataset('/work/opa/mi19918/Projects/nemo/BS/MOC_data/moc_meridional_sigma2_BSe3r1_mean.nc')
