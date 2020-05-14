from matplotlib import path


lon_region1 = np.array([ 0.7632,    8.1791,   14.5322,   21.6440,   26.3827,
                       24.2668, 17.4635,    8.0796,    4.0217,    0.7632]);
lat_region1 = np.array([81.0571,   82.7292,   83.5158,   83.5447,   82.5843,
                        81.2011, 80.8200,   79.2105,   79.1073,   81.0571]);


lons_lats_vect = np.column_stack((lon_region1, lat_region1)) # Reshape coordinates

p = path.Path(lons_lats_vect)

root_dir = '/tos-project1/NS9252K/CMIP6/ssp585/IPSL-CM6A-LR/r1i1p1f1/'
fname = root_dir + 'so_Omon_IPSL-CM6A-LR_ssp585_r1i1p1f1_gn_201501-210012.nc'

dfs = xr.open_dataset(fname)

root_dir = '/tos-project1/NS9252K/CMIP6/historical/IPSL-CM6A-LR/r1i1p1f1/'
fname = root_dir + 'so_Omon_IPSL-CM6A-LR_historical_r1i1p1f1_gn_195001-201412.nc'

dfc = xr.open_dataset(fname)

lon = np.copy(dfs.nav_lon)
lat = np.copy(dfs.nav_lat)
lons_lats_model = np.column_stack((lon.flatten(),lat.flatten()))
mask = p.contains_points(lons_lats_model)
mask.shape = lon.shape
mask = np.multiply(mask, 1)

# last 20 years of simulations 2081-2100
dfs = dfs.so[-240:,:,:,:]
# years between 1981-2000
dfc = dfc.so[(1980-1950+1)*12:(1980-1950+1)*12+240,:,:,:]
dfs = dfs.mean('time')
dfc = dfc.mean('time')

dfc.to_netcdf('IPSL_LR_salt_control.nc')
dfs.to_netcdf('IPLS_LR_salt_ssp585.nc')


for i in range(0,dfs.shape[2]):
    for j in range(0,dfs.shape[1]):
        if mask[j,i] == 1:
            plt.plot(dfs[:,j,i]-dfc[:,j,i],-dfs.olevel,'silver')






