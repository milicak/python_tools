from matplotlib import path


lon_region1 = np.array([ 0.7632,    8.1791,   14.5322,   21.6440,   26.3827,
                       24.2668, 17.4635,    8.0796,    4.0217,    0.7632]);
lat_region1 = np.array([81.0571,   82.7292,   83.5158,   83.5447,   82.5843,
                        81.2011, 80.8200,   79.2105,   79.1073,   81.0571]);


lons_lats_vect = np.column_stack((lon_region1, lat_region1)) # Reshape coordinates

p = path.Path(lons_lats_vect)

fname = '/tos-project1/NS9252K/CMIP6/ssp585/GFDL-CM4/r1i1p1f1/so_Omon_GFDL-CM4_ssp585_r1i1p1f1_gn_207501-209412.nc'
dfs = xr.open_dataset(fname, chunks={'time':10,'x': 100, 'y': 100})

fname = '/tos-project1/NS9252K/CMIP6/historical/GFDL-CM4/r1i1p1f1/so_Omon_GFDL-CM4_historical_r1i1p1f1_gn_199001-200912.nc'
dfc = xr.open_dataset(fname, chunks={'time':10,'x': 100, 'y': 100})

lon = np.copy(dfs.lon)
lat = np.copy(dfs.lat)
lons_lats_model = np.column_stack((lon.flatten(),lat.flatten()))
mask = p.contains_points(lons_lats_model)
mask.shape = lon.shape
mask = np.multiply(mask, 1)

# last 20 years of simulations 2081-2100
dfs = dfs.so[-240:,:,:,:]
# years between 1981-2000
dfc = dfc.so[(1980-1850+1)*12:(1980-1850+1)*12+240,:,:,:]
dfs = dfs.mean('time')
dfc = dfc.mean('time')

dfc.to_netcdf('GFDL_salt_control.nc')
dfs.to_netcdf('GFDL_salt_ssp585.nc')


for i in range(0,dfs.so.shape[2]):
    for j in range(0,dfs.so.shape[1]):
        if mask[j,i] == 1:
            plt.plot(dfs.so[:,j,i]-dfc.so[:,j,i],-dfs.lev,'r')
            # plt.plot(dfs[:,j,i]-dfc[:,j,i],-dfs.lev*1e-2,'silver')





