import numpy as np
import xarray as xr
import xesmf as xe


root_folder = '/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/OUT/'

gr = xr.open_dataset(root_folder + 'ocean_geometry.nc')

momNA12= xr.Dataset()
momNA12['lon'] = gr['geolon']
momNA12['lat'] = gr['geolat']

dft = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/WOA13/woa13_decav_t00_04v2.nc',decode_times=False)
dfs = xr.open_dataset('/Volumes/A1/workdir/milicak/datasets/WOA13/woa13_decav_s00_04v2.nc',decode_times=False)
dfwoa= xr.Dataset()
dfwoa['lon'] = dfs['lon']
dfwoa['lat'] = dfs['lat']


# Calculate remapping weights
# Using nearest neighbor - other options could be used here , e.g. bilinear.
regrid_woa = xe.Regridder(dfwoa, momNA12, 'patch',
                       periodic=False, 
                       filename='/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/regrid_woa_clim.nc',reuse_weights=True)


sst_woa = regrid_woa(dft['t_an'][0,0,:,:])   
sss_woa = regrid_woa(dfs['s_an'][0,0,:,:])   

sst_woa=sst_woa.to_dataset(name='sst_woa')
sss_woa=sss_woa.to_dataset(name='sss_woa')

sst_woa.to_netcdf('/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/sst_woa_NA12.nc')
sss_woa.to_netcdf('/Volumes/A1/workdir/milicak/datasets/MOM6/NA12/sss_woa_NA12.nc')


# ax.set_facecolor("gray")
# plt.pcolormesh(ssttemp-np.copy(sst_woa),cmap='RdBu_r',vmin=-4,vmax=4);plt.colorbar();

